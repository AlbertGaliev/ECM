using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace ECM2018.core
{
    interface IPoint
    {
        BigInteger x { get; }
        BigInteger y { get; }
        BigInteger z { get; }
        IPoint Multiply(BigInteger n);
        IPoint Add(IPoint P);
        IPoint Double();
        IPoint Clone();
        void ToPointAffine();
        bool IsInfinite();
        bool IsOnCurve();
    }
}
