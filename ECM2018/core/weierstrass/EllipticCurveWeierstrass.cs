using ECM2018.core.algorithms;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace ECM2018.core.weierstrass
{
    /// <summary>
    /// y^2=x^3+ax+b
    /// </summary>
    class EllipticCurveWeierstrass : IEllipticCurve
    {
        public BigInteger a { get; private set; }
        public BigInteger b { get; private set; }
        public BigInteger n { get; private set; }

        public EllipticCurveWeierstrass(BigInteger a, BigInteger b, BigInteger n)
        {
            this.a = a;
            this.b = b;
            this.n = n;
        }

        public static bool operator ==(EllipticCurveWeierstrass E1, EllipticCurveWeierstrass E2)
        {
            return E1.a == E2.a && E1.b == E2.b && E1.n == E2.n;
        }
        public static bool operator !=(EllipticCurveWeierstrass E1, EllipticCurveWeierstrass E2)
        {
            return !(E1 == E2);
        }
        public override bool Equals(Object obj)
        {
            if (obj == null)
                return false;
            EllipticCurveWeierstrass E = obj as EllipticCurveWeierstrass;
            return (new EllipticCurveWeierstrass(this.a, this.b, this.n) == E);
        }
        public override int GetHashCode()
        {
            return base.GetHashCode();
        }
        public static EllipticCurveWeierstrass GenerateCurveRandom(BigInteger n, out BigInteger x, out BigInteger y, out BigInteger g, bool aIsNull)
        {
            BigInteger a, b;
            do
            {
                do
                {
                    x = Algorithms.RandomBigInteger(n) % n;
                    y = Algorithms.RandomBigInteger(n) % n;
                    if (aIsNull)
                        a = 0;
                    else
                        a = Algorithms.RandomBigInteger(n) % n;
                    b = (y * y - x * x * x - a * x) % n;
                } while ((4 * a * a * a + 27 * b * b) % n == 0);
                g = Algorithms.GCD(n, 4 * a * a * a + 27 * b * b);
            } while (n == g);
            return new EllipticCurveWeierstrass(a, b, n);
        }


    }

}
