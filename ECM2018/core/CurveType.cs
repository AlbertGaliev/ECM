using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ECM2018.core
{
    enum TorsionType : byte
    {
        RandomCurve,
        Curve6,
        Curve12,
        Curve2x8,
        Curve12FromFile,
        Curve2x8FromFile
    }

    enum CurveType : byte
    {
        WeierstrassCurve,
        MontgomeryCurve,
        TwistedEdwardsCurve
    }

    enum CoordinatesType : byte
    {
        Affine,
        Projective,
        Jacobian,
        Inverted,
        Extended,
        Completed
    }
}
