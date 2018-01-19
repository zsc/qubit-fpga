module Deutsch where

import CLaSH.Prelude
import CLaSH.Sized.Fixed
import qualified Data.List as L

type RR = SFixed 1 10

--- Complex number utilities ---

type CC = Vec 2 RR

c0 = 0 :> 0 :> Nil
c1 = 1 :> 0 :> Nil
ci = 0 :> 1 :> Nil

sqr_norm :: CC -> RR
sqr_norm (a :> b :> Nil) = a * a + b * b

cadd :: CC -> CC -> CC
cadd = zipWith (+)

cneg :: CC -> CC
cneg (a :> b :> Nil) = -a :> -b :> Nil

csub :: CC -> CC -> CC
csub ca cb = cadd ca (cneg cb)

cmul :: CC -> CC -> CC
cmul (a :> b :> Nil) (c :> d :> Nil) = (a * c - b * d) :> (a * d + b * c) :> Nil

dotProduct xs ys = foldr cadd c0 (zipWith cmul xs ys)
matrixVector m v = map (`dotProduct` v) m

--- Qubit utilities ---

-- | amplitudes for |0> and |1>
type QBit = Vec 2 CC

q0 :: Signal QBit
q0 = register (c1 :> c0 :> Nil) q0

q1 :: Signal QBit
q1 = register (c0 :> c1 :> Nil) q1

qPlus = hadamardG q0
qMinus = hadamardG q1 

hadamard :: QBit -> QBit
hadamard = matrixVector ((h :> h :> Nil) :> (h :> (cneg h) :> Nil) :> Nil)
  where h = ($$(fLit (1 / sqrt 2)) :: RR) :> 0 :> Nil

hadamardG :: Signal QBit -> Signal QBit
hadamardG = register (repeat c0) . liftA hadamard

measure :: Signal QBit -> Signal RR
measure = register 0 . liftA (\ x ->  sqr_norm (x !! 1))

explode :: Signal QBit -> Signal QBit -> Signal (Vec 4 CC)
explode qx qy = register (repeat c0) $ liftA2 outer qx qy
  where 
    outer :: QBit -> QBit -> Vec 4 CC
    outer (x0 :> x1 :> Nil) y = (map (cmul x0) y) ++ (map (cmul x1) y)

measure0 :: Signal (Vec 4 CC) -> Signal RR
measure0 = register 0 . liftA (\ x -> sqr_norm (x !! 0) + sqr_norm (x !! 1))

make_complex = map (map (\ x -> x :> 0 :> Nil))

deutsch_u :: Vec 2 RR -> Vec 4 CC -> Vec 4 CC
deutsch_u (f0 :> f1 :> Nil) =
  matrixVector (make_complex (
                ((1 - f0) :> f1 :> 0 :> 0 :> Nil) :> 
                (f0 :> (1 - f1) :> 0 :> 0 :> Nil) :>
                (0 :> 0 :> (1 - f0) :> f1 :> Nil) :>
                (0 :> 0 :> f0 :> (1 - f1) :> Nil) :> Nil))

hadamard_I :: Vec 4 CC -> Vec 4 CC
hadamard_I =
  matrixVector (make_complex (
                (h :> 0 :> h :> 0 :> Nil) :> 
                (0 :> h :> 0 :> h :> Nil) :>
                (h :> 0 :> - h :> 0 :> Nil) :>
                (0 :> h :> 0 :> - h :> Nil) :> Nil))
  where h = $$(fLit (1 / sqrt 2)) :: RR

deutsch :: Vec 2 RR -> Signal RR
deutsch f0f1 =
  let xy = explode qPlus qMinus in
  let xy2 = register (repeat c0) $ liftA (deutsch_u f0f1) xy in
  let xy3 = register (repeat c0) $ liftA hadamard_I xy2 in
  measure0 xy3

topEntity :: Signal (Vec 4 RR)
topEntity = bundle (map deutsch (f0 :> f1 :> f2 :> f3 :> Nil))
  where f0 = 0 :> 0 :> Nil
        f1 = 1 :> 1 :> Nil
        f2 = 0 :> 1 :> Nil
        f3 = 1 :> 0 :> Nil
