(ns codon-bias.cai
  (:require [clojure.core.reducers :as r])
  (:use edu.bc.bio.seq-utils
        [edu.bc.utils :only [log]]
        [edu.bc.utils.probs-stats :only [mean]]))

;;;-------For testing purposes

(def testseqs ["TTCACGTCTTGCTCTGCGACCTCTTTAACG" "GGGGATAACGTGCCGCTATGTGGCACGCAA"
               "AATTATCCGTTTGACCCCAAGGCGCGTATT" "TCCCATCCGTCCTGTGCATCTGCAGAACCA"
               "CGGTTACGAGTAGCGCATAGCCCCTGCCGA" "ACAGTTGTACTCGCCCTCCCCCATAAAACG"
               "GTTAGAAGGAGAGTTAGTCCTAAACAGTGG" "GCGGCGAGCGAACGCCGTCACACAACGCGC"
               "ACCCTGGTCCAAACCAAACATTAAGCCGAA" "ACTCCGGGCAAGTCGCTGGTCGCCAGAACA"])

(def test-background-seqs
  (repeatedly 100000 #(->> (interleave [\A \C \G \T] (repeat 0.25))
                       (apply hash-map )
                       (generate-rand-seq 30))))

(def testseq "GAGGTTCAAGAATTGGGGCATAGAACGGAC")
;;;-------------------------

;;;---utility functions
(defn str-partition-all [n s] (map #(apply  str %) (partition-all n s)))

(defn remove-start-stop [sq]
  (second (re-find #"AUG(.*)(TGA|TAG|TAA)" sq)))
;;;-----------------------------

(defn get-seqs [query])

(defn process-seqs
  "Reads in sequences to be used for a background set. Returns a map
  of amino acids, codons and respective counts"
  
  [inseqs]
  (->> (map #(frequencies (str-partition-all 3 %)) inseqs)
       (apply merge-with +)
       (reduce (fn [M [codon cnt]]
                 (let [aa (+DNA-CONDON-MAP+ codon)
                       V (get M aa [])]
                   (assoc M aa (conj V [codon cnt]))))
               {})))

(defn calc-rscu
  "Calculates the rscu for synonymous codons of a given amino acid."

  [aa syn-codons]
  (let [d (+REVERSE-DNA-CODON-MAP+ aa)
        syn-codons (merge (apply hash-map (interleave d (repeat 0.5)));add pseudocounts
                          (into {} syn-codons))]
    (reduce (fn [M codon]
              (assoc M codon (/ (get syn-codons codon 0.5)
                                (mean (vals syn-codons)))))
            {} d)))

(defn calc-weights
  "Calculates the weights of codons given a set of sequences"

  [sqs]
  (let [aa-counts (-> (process-seqs sqs)
                      (dissoc "" "PY" "SE"))] ;remove stop codons
    (reduce 
     (fn [M [aa codons]]
       (let [rscu (calc-rscu aa codons)
             max-rscu (apply max (vals rscu))]
         (->> rscu
              (map #(vector (first %) (/ (second %) max-rscu)))
              (into {})
              (merge M))))
     {}  aa-counts)))

(defn calc-cai
  "Given a sequence and weights of background sequences, calculates the CAI for a given sequence"

  [weights sq]
  (let [sq (str-partition-all 3 sq)]
    (->> (map weights sq)
         (remove nil?)
         (map log)
         mean
         (Math/exp))))

;;;-----------tests
(defn test-boundary []
  (let [sq (fn [] (generate-rand-seq 30 (->> (interleave [\A \C \G \T] (repeat 0.25))
                                           (apply hash-map ))))]
    (zero?
     (->> (for [s (repeatedly 10000 sq)]
            (calc-cai (calc-weights test-background-seqs) s))
          (remove #(< 0 % 1))
          count))))

(defn test-cai []
  (prn "should be approx -0.6173786103910937")
  (calc-cai (calc-weights testseqs)
            (apply str '("GAG" "GTT" "CAA" "GAA" "TTG" "GGG" "CAT" "AGA" "ACG" "GAC")))
  )
