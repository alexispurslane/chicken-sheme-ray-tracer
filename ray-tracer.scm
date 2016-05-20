(declare (standard-bindings)
         (extended-bindings)
         (block)
         (not safe))
(require-extension list-of)
(define pi 3.14159265358979323846264338328)
(define-record vec3 x y z)
(define (vec3f x) (vec3 x x x))

(define-syntax def-vec3-op
  (er-macro-transformer
   (lambda (exp rename compare)
     (let ([name (cadr exp)] [acc  (caddr exp)] [op (cadddr exp)])
       `(define (,(string->symbol (string-append
                                   "vec3:" (symbol->string name))) x y)
          (,acc (,op (vec3-x x) (vec3-x y))
                (,op (vec3-y x) (vec3-y y))
                (,op (vec3-z x) (vec3-z y))))))))

(def-vec3-op dot + *)
(def-vec3-op - vec3 -)
(def-vec3-op + vec3 +)
(def-vec3-op * vec3 *)
(def-vec3-op / vec3 /)
(define (vec3:dot v) (vec3 (- (vec3-x v)) (- (vec3-y v)) (- (vec3-z v))))

(def-vec3-op < and <)
(def-vec3-op > and >)
(def-vec3-op >= and >=)
(def-vec3-op <= and <=)
(def-vec3-op = and =)

(define (vec3:dotength2 v)
  (+ (expt (vec3-x v) 2) (expt (vec3-y v) 2) (expt (vec3-z v) 2)))

(define (vec3:dotormalize v)
  (let ([nor2 (vec3:dotength2 v)]
        [invNor (/ 1 (sqrt nor2))])
    (if (not (zero? nor2))
        (vec3 (* (vec3-x v) invNor)
              (* (vec3-y v) invNor)
              (* (vec3-z v) invNor))
        v)))

(define-record ray origin dir)

(define-record sphere
  center radius radius2
  surface-color emission-color
  transparency reflection)

(define (intersect s r)
  (let* ([l   (- (sphere-center s) (ray-origin r))]
         [tca (vec3:dot l (ray-dir r))]
         [d2  (- (vec3:dot l l) (* tca tca))]
         [thc (sqrt (- radius2 d2))])
    (if (or (negative? tca)
            (> d2 radius2))
        (values #f 0 0)
        (values #t (- tca thc)
                   (+ tca thc)))))

(define (mix a b mix) (+ (* b mix) (* a (- 1 mix))))

(define MAX-RAY-DEPTH 5)

(define (trace view-ray spheres depth)
  (let-values ([(tnear sphere) (foldl (lambda (nearest sphere)
                                        (let-values ([(intersecting? t0 t1) (intersect sphere view-ray)])
                                          (if intersecting?
                                              (let-values ([ct (if (negative? t0) t1 t0)]
                                                           [(tnear cnear) nearest])
                                                (if (< ct tnear)
                                                    (values ct sphere)
                                                    (values tnear cnear))))))
                                      (values +inf.0 (void)) spheres)])
    (if (void? sphere)
        (vec3f 2.0)
        (let* ([phit (vec3:dot (ray-origin view-ray)
                               (vec3:dot (ray-dir view-ray)
                                         (vec3 tnear tnear tnear)))]
               [nhit (vec3:dotormalize (vec3:dot phit (sphere-center sphere)))]
               [bias 1e-4]
               [inside? (positive? (vec3:dot (ray-dir view-ray) nhit))]
               [nhit (if inside? (vec3:dot nhit) nhit)])
          (if (and (or (positive? (sphere-transparency sphere))
                       (positive? (sphere-reflection sphere)))
                   (< depth MAX-RAY-DEPTH))
              (let* ([facing-ratio (vec3:dot (vec3:dot (ray-dir view-ray) nhit))]
                     [fresnel-effect (mix (expt (- 1 facing-ratio) 3) 1 0.1)]
                     [refl-dir (vec3:dotormalize
                                (* (vec3:dot (ray-dir ray-view) nhit) 2
                                   (vec3:dot (ray-dir ray-view) nhit)))]
                     [reflection (trace
                                  (make-ray (vec3:dot phit (* nhit bias)) refl-dir)
                                  spheres
                                  (+ depth 1))]
                     [refraction (if (sphere-transparency sphere)
                                     (let* ([ior 1.1]
                                            [eta (if inside? ior (/ 1 ior))]
                                            [cosi (vec3:dot
                                                   (vec3:dot nhit
                                                             (ray-dir view-ray)))]
                                            [k (- 1 (* eta eta (- 1 (* cosi cosi))))]
                                            [refr-dir (vec3:dotormalize
                                                       (+ (* (ray-dir view-ray) eta)
                                                          (* nhit eta (- cosi (sqrt k)))))])
                                       (trace (make-ray (vec3:dot phit (* nhit bias))
                                                   refr-dir)
                                              spheres (+ depth 1)))
                                     (vec3f 0.0))])
                (vec3:dot (vec3:dot reflection fresnel-effect)
                          (vec3:dot refraction (- 1 fresnel-effect)
                                    (sphere-transparency sphere)
                                    (sphere-surface-color sphere))))
              (+ (sphere-emission-color sphere)
                 (foldl (lambda (surface-color spherei)
                          (if (> (vec3-x (sphere-emission-color sphere)) x)
                              (let* ([light-direction (normalize
                                                       (- (sphere-center spherei)
                                                          phit))]
                                     [transmission
                                      (foldl
                                       (lambda (ts spherej)
                                         (if (not (equal? spherej spherei))
                                             (let-values ([(i? t0 t1) (intersect
                                                                       spherej
                                                                       (make-ray
                                                                        (+ phit (* nhit bias))
                                                                        light-direction))]
                                                          [(t stop?) ts])
                                               (if (or i? stop?)
                                                   (values 0 #t)
                                                   (values t stop?)))))
                                       1 spheres)])
                                (+ surface-color (* (sphere-surface-color sphere)
                                                    transmission
                                                    (max 0 (vec3:dot nhit light-direction))
                                                    (sphere-emission-color spherei))))))
                        0 spheres)))))))

(define (real->u8 n)
  (inexact->exact
   (max 0 (min 255 (if (integer? n) n (floor n))))))

(define write-u8 write-byte)

(define (write-color v1)
  (write-u8 (real->u8 (* (vec3-x v1) 255)))
  (write-u8 (real->u8 (* (vec3-y v1) 255)))
  (write-u8 (real->u8 (* (vec3-z v1) 255))))

(define (render spheres #!optional [width 640] [height 480])
  (let* ([inv-width (/ 1 width)] [inv-height (/ 1 height)]
         [fov 30] [aspect-ratio (/ width height)]
         [angle (tan (/ (* pi 0.5 fov) 180))]
         [image (map (lambda (xs)
                       (map (lambda (xy)
                              (let* ([x (car x)] [y (cadr y)]
                                     [xx (* (- (* 2 (+ x 0.5 ) inv-width) 1)
                                            angle
                                            aspect-ratio)]
                                     [yy (* (- 1 (* 2 (+ y 0.5) inv-height)) angle)]
                                     [ray-dir (vec3:normalize (vec3 xx yy -1))])
                                (trace (make-ray (vec3f 0) ray-dir) spheres 0))) xs))
                     (list-of (list-of (list x y) (x range 0 width)) (y range 0 height)))])
    (with-output-to-file "image.ppm"
      (lambda ()
        (print "Writing...")
        (display-args "P6\n")
        (display-args width " " height "\n")
        (display-args "255\n")
        (for-each (lambda (c) (write-color c)) image)))))
