using GroveAlg
using Base.Test

d = 5
@test TreeCheck(3,TreeIndex([1,2,3])) && GroveCheck([1,2,3])
@test sum(GroveError(d)) == 0 == sum(GroveError(Grove(d))) == sum(GroveError(Grove(d).Y))
@test TreeCheck([1,2,3]) == TreeCheck(PBTree([1,2,3])) == TreeCheck(Grove(3)) == TreeCheck(Grove(3).Y) == PBTree(3,5) |> TreeCheck
@test TreeIndexCn(d) == [i//Cn(d) for i ∈ 1:Cn(d)] == TreeIndex(Grove(d).Y)//Cn(d)
@test σ(σ(Grove(3,7))) == Grove(3,7)
@test Graft([1,2,3],[2,1]) |> Grove == Graft([1,2],[3,2,1]) |> σ
@test GroveAlg.PrimitiveTree([1,2,3]) && GroveAlg.PrimitiveTree([3,2,1])
@test GroveAlg.BranchLeft([1,4,2,1]) == GroveAlg.BranchRight([1,2,4,1])
@test GroveAlg.LeftInherited([1,2,3]) == GroveAlg.RightInherited([3,2,1])
@test [1,2,3] |> GroveSort! == [1,2,3] |> Grove
@test Grove(d) |> GroveSort! == Grove(d) == Grove(d,2^Cn(d)-1)
@test Grove(d) == d |> GroveAlg.TreeBase |> GroveAlg.TreeLoday
@test (k = Grove(d) |> GroveAlg.TreeBase |> GroveAlg.TreeLoday;
  k == k |> GroveAlg.TreeBase |> GroveAlg.TreeLoday)
@test [GroveAlg.TreeBase([1:d...])] == GroveAlg.TreeBase(d,1)
@test Grove(d) |> GroveAlg.TreeBase |> TreeInteger == Grove(d) |> TreeInteger
@test (GroveSort(false); GroveAlg.TreeLoday(d,TreeIndex(Grove(d))) == Grove(d))
@test (GroveSort(false); Grove(d)) == (GroveSort(true); Grove(d))
@test Grove(d,GroveIndex(Grove(d).Y)) == GroveAlg.TreeLoday(d)
@test GroveBin([1,4,2,1]) |> Grove |> GroveCheck
@test GroveAlg.TreeLoday(d,GroveBit(Grove(d))) == Grove(d) |> GroveBin |> Grove
@test Grove(d) |> GroveIndex == 2^Cn(d)-1 == Grove(d) |> GroveBit |> GroveIndex
@test TreeRational(d) == TreeRational(Grove(d))
@test TreeRational(d,BigInt(d)) == Grove(d,d) |> TreeRational
@test (TreeShift(false); TreeRational(d,GroveBit(Grove(d))) == TreeRational(Grove(d)))
@test (TreeShift(true); TreeRational(d,GroveBit(Grove(d).Y)) == TreeRational(Grove(d)))
@test TreeInteger(d,d)==TreeInteger(d,GroveBit(d))==TreeInteger(d,TreeIndex(Grove(d,d)))
@test (j = ([1,2]⊣[1,2])+[1,2,3]; GroveAlg.TreeLoday(7,GroveIndex(j)) == j)
@test (j = [1,2,3]+([1,2]⊢[1,2]); GroveAlg.TreeLoday(7,GroveIndex(j)) == j)
@test (Grove([1,2,3])⊣PBTree([1,3,1])) == (PBTree([1,2,3])⊣Grove([1,3,1])) == (Grove([1,2,3])⊣Grove([1,3,1])) == ((Grove([1,2,3]).Y)⊣(PBTree([1,3,1]).Y))
@test (Grove([1,2,3])⊢PBTree([1,3,1])) == (PBTree([1,2,3])⊢Grove([1,3,1])) == (Grove([1,2,3])⊢Grove([1,3,1])) == ((Grove([1,2,3]).Y)⊢(PBTree([1,3,1]).Y))
@test Grove(8,GroveIndex(Grove(5,1000)+Grove(3,7)))==Grove(5,1000)+Grove(3,7)
@test GroveComposition(3) == GroveComposition(3)
@test GrovePrint(Grove(3,7)) == GrovePrint(3) == GrovePrint(GroveAlg.TreeBase([1:d...]))
@test GroveAlg.TreeBase(d,1) |> GrovePrint == GrovePrint([Grove(d,1)])
