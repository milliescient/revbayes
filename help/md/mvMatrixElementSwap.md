## name
mvMatrixElementSwap
## title
Matrix Element Swap Move
## description
This move proposes exchanging two elements of a matrix.

## details
The move picks two elements at random and exchanges them. By default, any two elements can be exchanged.

The `margin` argument follows R's MARGIN convention: 1 swaps within a random row, and 2 within a random column. You can instead specify a row with `row` or a column with `col`, counting from 1. Give at most one of `margin`, `row` and `col`.
## authors
June Walker
## see_also
mvMatrixElementScale
mvMatrixElementSlide
## example
    moves = VectorMoves()
    x ~ dnDecomposedInvWishart(df=10, kappa=2.0, dim=4)
    # any two elements of the matrix
    moves.append( mvMatrixElementSwap(x, weight=1) )
    # two elements of column 2
    moves.append( mvMatrixElementSwap(x, col=2, weight=1) )

## references
