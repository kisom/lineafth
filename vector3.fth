\ Copyright (c) 2016 Kyle Isom <coder@kyleisom.net>
\ 
\ Permission is hereby granted, free of charge, to any person obtaining a
\ copy of this software and associated documentation files (the "Software"),
\ to deal in the Software without restriction, including without limitation
\ the rights to use, copy, modify, merge, publish, distribute, sublicense,
\ and/or sell copies of the Software, and to permit persons to whom the
\ Software is furnished to do so, subject to the following conditions:
\ 
\ The above copyright notice and this permission notice shall be included
\ in all copies or substantial portions of the Software.
\ 
\ THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
\ IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
\ FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
\ THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
\ OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
\ ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
\ OTHER DEALINGS IN THE SOFTWARE.

\ vector3.fth: definitions for R^3 vectors.


( Compute the integer square root of a number. Note that this computes )
( the floor of the square root. For example, 14 ISQRT => 3, but        )
( 14 SQRT = 3.74.                                                      )
( n -- n )
: ISQRT   -1 BEGIN 1+ 2DUP DUP * - 0< UNTIL 1- ;

( Square and cubing of numbers.                                        )
( n -- n )
: ^2   DUP * ;
: ^3   DUP ^2 * ;


\ Various 3-tuple variants of stack manipulation functions.

( w1 w2 w3 -- w1 w2 w3 w1 w2 w3 )
: 3DUP   DUP 2OVER ROT ;
( w1 w2 w3 w4 w5 w6 -- w4 w5 w6 w1 w2 w3 )
: 3SWAP  5 ROLL 5 ROLL 5 ROLL ;
( w1 w2 w3 -- )
: 3DROP   2DROP DROP ;


( Add two vectors together.                                            )
( w1 w2 w3 w4 w5 w6 -- w7 w8 w9 )
: VEC3ADD
    3 ROLL +              ( add k components )
    SWAP 3 ROLL +         ( add j components )
    2SWAP +               ( add i components )
    SWAP ROT              ( restore order    )
;

( Compute the magnitude of a vector.                                   )
( w1 w2 w3 -- w4 )
: VEC3MAG   >R >R ^2 R> ^2 + R> ^2 + ISQRT SWAP DROP ;

( Scalar multiply the R3 vector. Note that the vector comes first,     )
( followed by the scalar. For example,                                 )
( 1 2 3 10 VEC3SMUL => 10 20 30                                        )
( This version does not use the return stack. I don't know if this is  )
( a useful feature or not.                                             )
( w1 w2 w3 w4 -- w5 w6 w7)
: VEC3SMUL   DUP ROT * SWAP 2SWAP ROT DUP ROT * SWAP ROT * SWAP ROT ;

( Scalar multiply the R3 vector. Note that the vector comes first,     )
( followed by the scalar. For example,                                 )
( 1 2 3 10 VEC3SMUL => 10 20 30                                        )
( w1 w2 w3 w4 -- w5 w6 w7)
( This version does use the return stack.                              )
: VEC3SMULR
    SWAP 2SWAP ROT    ( Push k below vectors.            )
    >R >R >R          ( Push components to return stack. )
    DUP R> * SWAP     ( Scale i component.               )
    DUP R> * SWAP     ( Scale j component.               )
    R> *              ( Scale k component.               )
;

( Pretty-print a vector.                                               )
( w1 w2 w3 -- w1 w2 w3 )
: SHOVEC3
    ." Vector: "
    >R >R DUP DUP 1 <> IF . ELSE DROP THEN ." i + "
    R>    DUP DUP 1 <> IF . ELSE DROP THEN ." j + " 
    R>    DUP DUP 1 <> IF . ELSE DROP THEN ." k " 
;


\ The following demo adds two vectors:
\     a = (1,2,3)
\     b = (4,5,6)
\ which should produce the resulting vector c = (5,7,9).
\ 1 2 3   CR SHOVEC3    ( vector a )
\ 4 5 6   CR SHOVEC3    ( vector b )
\ VEC3ADD CR SHOVEC3    ( vector c )


( Compute the dot product of two vectors.                              )
( w1 w2 w3 w4 w5 w6 -- w7 w8 w9 )
: VEC3DOT
    2SWAP SWAP ROT + >R   ( sum the k components                     )
    SWAP ROT + >R         ( sum the j componennts                    )
    + R> R>               ( sum the i components and restore j and k )
;
    
( Compute the cross-product of two vectors.                            )
( w1 w2 w3 w4 w5 w6 -- w7 w8 w9 )
: VEC3CROSS
    ( The comments use the following matrix:     )
    (     | i j k |                              )
    (     | a b c |                              )
    (     | x y z |                              )
    ( The vector <a, b, c> is the first          )
    ( vector on the stack, and <x, y, z> is      )
    ( the second vector. E.g.                    )
    ( a b c x y z -- i j k                       )
    ( [bz - cy]i + [cx - az]j + [ay - bx]k       )

    ( Compute the i component.                   )
    DUP 5 ROLL DUP ROT * >R          ( b*z       )
    4 ROLL DUP 4 ROLL DUP ROT *      ( c*y       )
    R> SWAP - >R                     ( b*z - c*y )

    ( Compute the j component.                   )
    4 ROLL DUP 3 ROLL * >R           ( c*x       )
    4 ROLL DUP 5 ROLL *              ( a*z       )
    R> SWAP - >R                     ( c*z - a*z )

    ( Compute the k component.                   )
    ROT * ROT ROT * -                ( a*y - b*x )

    ( Restore the i and j components, then       )
    ( reverse the resulting vector.              )
    R> R> SWAP ROT
;

\ Cross-product test vectors:
\ The first comes from the Octave manual, which was used to verify
\ both answers.
\ 1 1 0 0 1 1 VEC3CROSS CR SHOVEC3 CR
\  should produce (1, -1, 1)
\ 1 2 3 4 5 6 VEC3CROSS CR SHOVEC3 CR
\  should produce (-3, 6, -3)

( Roll a vector on top of the stack. This is used to push a scale      )
( value down from the stack.                                           )

: VEC3ROLL   3 ROLL 3 ROLL 3 ROLL ;    

( Compute the unit vector of an R3 vector.                             )
( w1 w2 w3 w4 -- w5 w6 w7 )
\ note: in progress
\ : VEC3UNIT
\     VEC3ROLL                   ( push scaler off stack )
\     3DUP VEC3MAG               ( compute magnitude     )
\     DUP -ROT 5 ROLL DUP VEC3ROLL -ROT */ >R
\                                ( scale divide k by     )
\                                (   magnitude           )
\     SWAP DUP 3 ROLL 3 ROLL DUP VEC3ROLL ROT */ >R
\                                ( scale divide j by     )
\                                (   magnitude           )
\     -ROT */                    ( scale divide i by     )
\                                  (   magnitude           )
\     R> R>                      ( restore i and j       )
\ ;
