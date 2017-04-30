using GroveAlg
using Base.Test

# write your own tests here
d = 5
@test TreeCheck(3,TreeIndex([1,2,3])) && GroveCheck([1,2,3])
@test sum(GroveError(d)) == 0
@test σ(σ(Grove(3,7))) == Grove(3,7)
@test Graft([1,2,3],[2,1]) == σ(Graft([1,2],[3,2,1]))
@test GroveAlg.PrimitiveTree([1,2,3]) && GroveAlg.PrimitiveTree([3,2,1])
@test GroveSort!([1,2,3]) == Grove([1,2,3])
@test GroveSort!(Grove(d)) == Grove(d) == Grove(d,2^Cn(d)-1)
@test Grove(d) == GroveAlg.TreeLoday(GroveAlg.TreeBase(Grove(d)))
@test (k = GroveAlg.TreeLoday(GroveAlg.TreeBase(Grove(d))); k == GroveAlg.TreeLoday(GroveAlg.TreeBase(k)))
@test TreeInteger(GroveAlg.TreeBase(Grove(d))) == TreeInteger(Grove(d))
@test (GroveSort(false); GroveAlg.TreeLoday(d,TreeIndex(Grove(d))) == Grove(d))
@test (GroveSort(false); Grove(d)) == (GroveSort(true); Grove(d))
@test Grove(d,GroveIndex(Grove(d))) == Grove(d)
@test GroveAlg.TreeLoday(3,GroveBit(Grove(3))) == Grove(GroveBin(Grove(3)))
@test GroveIndex(Grove(d)) == 2^Cn(d)-1
@test (TreeShift(false); TreeRational(d,GroveBit(Grove(d))) == TreeRational(Grove(d)))
@test (TreeShift(true); TreeRational(d,GroveBit(Grove(d))) == TreeRational(Grove(d)))
@test (j = ([1,2]⊣[1,2])+[1,2,3]; Grove(7,GroveIndex(j)) == j)
@test Grove(8,GroveIndex(Grove(5,1000)+Grove(3,7)))==Grove(5,1000)+Grove(3,7)
@test GroveComposition(3) == GroveComposition(3)
@test GrovePrint(Grove(3,7)) == GrovePrint(Grove(3,4))
