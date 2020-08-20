/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.6                                   |
|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    `format'      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// block definition for a porosity with an angled inlet/outlet
// the porosity is not aligned with the main axes
//
dnl> -----------------------------------------------------------------
dnl> <STANDARD DEFINTIONS>
dnl>
changecom(//)changequote([,]) dnl>
define(calc, [esyscmd(perl -e 'print ($1)')]) dnl>
define(VCOUNT, 0)  dnl>
define(vlabel, [[// ]pt VCOUNT ($1) define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])  dnl>
dnl>
define(hex2D, hex ($1b $2b $3b $4b $1f $2f $3f $4f)) dnl>

define(quad2D, ($1f $1b $2b $2f))  dnl>

define(frontQuad, ($1f $2f $3f $4f)) dnl>
define(backQuad, ($4b $3b $2b $1b)) dnl>
dnl>
dnl> </STANDARD DEFINTIONS>
dnl> -----------------------------------------------------------------
dnl>

//refine times
define(refine,3)   

define(n1, calc(refine*30)) dnl>
define(n2, calc(refine*50)) dnl>
define(n3, calc(refine*100)) dnl>
define(n4, calc(refine*50)) dnl>
define(n5, calc(refine*10)) dnl>
define(n6, calc(refine*12)) dnl>

define(L1,0.3) //domain length 1
define(L2,0.5) //domain length 2
define(L3,1) //domain length 3
define(L4,0.405) //domain length 3
define(L5,0.025) //domain length 3
define(L6,0.12) //domain length 3

define(x0,calc(L1)) dnl>
define(x1,calc(L1+L2)) dnl>
define(x2,calc(L1+L2+L3)) dnl>
define(x3,0) dnl>
define(x4,calc(L1)) dnl>
define(x5,calc(L1+L2)) dnl>
define(x6,calc(L1+L2+L3)) dnl>
define(x7,0) dnl>
define(x8,calc(L1)) dnl>
define(x9,calc(L1+L2)) dnl>
define(x10,calc(L1+L2+L3)) dnl>

define(y0,0) dnl>
define(y1,0) dnl>
define(y2,0) dnl>
define(y3,calc(L4)) dnl>
define(y4,calc(L4)) dnl>
define(y5,calc(L4)) dnl>
define(y6,calc(L4)) dnl>
define(y7,calc(L4+L5)) dnl>
define(y8,calc(L4+L5)) dnl>
define(y9,calc(L4+L5)) dnl>
define(y10,calc(L4+L5)) dnl>

define(zFront,calc(L6/2.0)) dnl>
define(zBack,calc(-L6/2.0)) dnl>

define(gp1,0.1) dnl>
define(gp2,5) dnl>

convertToMeters 1;

vertices
(
    //vertex on the front
    (x0   y0 zFront)  vlabel(a0f)
    (x1   y1 zFront)  vlabel(a1f)
    (x2   y2 zFront)  vlabel(a2f)
    (x3   y3 zFront)  vlabel(a3f)
    (x4   y4 zFront)  vlabel(a4f)
    (x5   y5 zFront)  vlabel(a5f)
    (x6   y6 zFront)  vlabel(a6f)
    (x7   y7 zFront)  vlabel(a7f)
    (x8   y8 zFront)  vlabel(a8f)
    (x9   y9 zFront)  vlabel(a9f)
    (x10  y10 zFront)  vlabel(a10f)

    //vertex on the back
    (x0   y0 zBack)  vlabel(a0b)
    (x1   y1 zBack)  vlabel(a1b)
    (x2   y2 zBack)  vlabel(a2b)
    (x3   y3 zBack)  vlabel(a3b)
    (x4   y4 zBack)  vlabel(a4b)
    (x5   y5 zBack)  vlabel(a5b)
    (x6   y6 zBack)  vlabel(a6b)
    (x7   y7 zBack)  vlabel(a7b)
    (x8   y8 zBack)  vlabel(a8b)
    (x9   y9 zBack)  vlabel(a9b)
    (x10  y10 zBack)  vlabel(a10b)
);

blocks
(
    // block 1
    hex2D(a0, a1, a5, a4)
    block1 ( n2 n4 n6 )  simpleGrading (1 gp1 1)

    // block 2
    hex2D(a1, a2, a6, a5)
    block1 ( n3 n4 n6 )  simpleGrading (1 gp1 1)

    // block 3
    hex2D(a3, a4, a8, a7)
    block1 ( n1 n5 n6 )  simpleGrading (1 gp2 1)

    // block 4
    hex2D(a4, a5, a9, a8)
    block1 ( n2 n5 n6 )  simpleGrading (1 gp2 1)

    // block 5
    hex2D(a5, a6, a10, a9)
    block1 ( n3 n5 n6 )  simpleGrading (1 gp2 1)

);

edges
(
);

patches
(
    symmetryPlane top
    (
       (a7f a8f a8b a7b)
       (a8f a9f a9b a8b)
       (a9f a10f a10b a9b)
    )

    patch inlet
    (
       (a3f a7f a7b a3b)
    )

    wall walls
    (
       (a3f a4f a4b a3b)
       (a0f a4f a4b a0b)
       (a0f a1f a1b a0b)
       (a2f a6f a6b a2b)
       (a10f a6f a6b a10b)
    )

    patch bottomInlet
    (
       (a1f a2f a2b a1b)
    )

    wall frontAndBack
    (
       (a0b a1b a5b a4b)
       (a0f a1f a5f a4f)

       (a1b a2b a6b a5b)
       (a1f a2f a6f a5f)

       (a3b a4b a8b a7b)
       (a3f a4f a8f a7f)

       (a4b a5b a9b a8b)
       (a4f a5f a9f a8f)

       (a5b a6b a10b a9b)
       (a5f a6f a10f a9f)
    )

);

mergePatchPairs
(
);

// ************************************************************************* //
