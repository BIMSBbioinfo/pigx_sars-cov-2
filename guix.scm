;;; PiGx SARS-CoV-2 wastewater sequencing pipeline
;;; Copyright © 2021, 2022 Ricardo Wurmus <rekado@elephly.net>
;;;
;;; This file is part of PiGx SARS-CoV-2 wastewater sequencing pipeline
;;;
;;; This is free software; see LICENSE file for details.

(use-modules
 (guix build-system gnu)
 (guix download)
 (guix packages)
 (guix licenses)
 (gnu packages)
 (ice-9 match))

(include "manifest.scm")

(define (p name)
  `(,name ,(specification->package name)))

(define %version
  (symbol->string (with-input-from-file "VERSION" read)))

(define (bimsb-origin name hash)
  (origin
    (method url-fetch)
    (uri
     (string-append "https://bimsbstatic.mdc-berlin.de/akalin/AAkalin_pathogenomics"
                    "/databases_small-20221006/" name))
    (sha256 (base32 hash))))

(define-public pigx-sars-cov-2-development
  (package
    (name "pigx-sars-cov-2")
    (version %version)
    (source
     (string-append (getcwd) "/pigx_sars-cov-2-" version ".tar.gz"))
    (build-system gnu-build-system)
    (arguments
     (list
      #:phases
      '(modify-phases %standard-phases
         (add-after 'unpack 'unpack-databases
           (lambda* (#:key inputs #:allow-other-keys)
             ;; The tests need to be able to write caches to HOME.
             ;; They also default to reading the databases from there.
             (setenv "HOME" "/tmp")
             ;; Unpack the three databases in the expected location.
             (let ((root "/tmp/.local/share/pigx/databases")
                   (use-underscore (lambda (c) (if (equal? c #\-) #\_ c))))
               (for-each (lambda (db)
                           (let ((where (string-append root "/"
                                                       (string-map use-underscore db))))
                             (mkdir-p where)
                             (invoke "tar" "-C" where
                                     "-xf" (assoc-ref inputs db))))
                         '("kraken-db" "krona-db" "vep-db"))))))))
    (inputs
     (map p %packages))
    (native-inputs
     `(("kraken-db"
        ,(bimsb-origin
          "kraken_db.tar.gz"
          "0sdm4xh5npg6c3y2pz8xgphim4qpglm8wdid6rlaaqsn6iikv0mz"))
       ("krona-db"
        ,(bimsb-origin
          "krona_db.tar.gz"
          "1rwy4gd3vw1gdjldrgf44c1xaa3vq8i3pgisjhrac81yx63x8f2h"))
       ("vep-db"
        ,(bimsb-origin
          "vep_db.tar.gz"
          "0d8hhi43zsw3wqm7gd0z0gpcdsc6h6ra0imn87hifl9a64jxqzxz"))
       ,@(map p %native-packages)))
    (home-page "https://bioinformatics.mdc-berlin.de/pigx/")
    (synopsis "SARS-CoV-2 wastewater sequencing pipeline")
    (description "")
    (license gpl3+)))

pigx-sars-cov-2-development
