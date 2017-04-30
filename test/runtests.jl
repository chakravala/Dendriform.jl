using GroveAlg
using Base.Test

# write your own tests here
d = 5
@test GroveSort!(Grove(d)) == Grove(d) == Grove(d,2^Cn(d)-1)
@test Grove(d) == GroveAlg.TreeLoday(GroveAlg.TreeBase(Grove(d)))
@test (k = GroveAlg.TreeLoday(GroveAlg.TreeBase(Grove(d))); k == GroveAlg.TreeLoday(GroveAlg.TreeBase(k)))
@test TreeInteger(GroveAlg.TreeBase(Grove(d))) == TreeInteger(Grove(d))
@test Grove(d,GroveIndex(Grove(d))) == Grove(d)
@test (j = ([1,2]⊣[1,2])+[1,2,3]; Grove(7,GroveIndex(j)) == j)
@test Grove(8,GroveIndex(Grove(5,1000)+Grove(3,7)))==Grove(5,1000)+Grove(3,7)
@test (GroveSort(false); Grove(d)) == (GroveSort(true); Grove(d))
