using ECM2018.core.algorithms;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace ECM2018.core.montgomery
{
    /// <summary>
    /// b(y^2)=x^3+a(x^2)+x
    /// </summary>
    class EllipticCurveMontgomery : IEllipticCurve
    {
        public BigInteger a { get; private set; }
        public BigInteger b { get; private set; }
        public BigInteger n { get; private set; }

        public BigInteger An { get; private set; }
        public BigInteger Ad { get; private set; }
        public BigInteger a24 { get; private set; }

        public EllipticCurveMontgomery(BigInteger An, BigInteger Ad, BigInteger n)
        {
            this.An = An;
            this.Ad = Ad;
            this.a = An / Ad - 2;
            this.n = n;
        }

        public EllipticCurveMontgomery(BigInteger a, BigInteger n)
        {
            this.a = a;
            this.b = (a + 2) / 4;
            this.n = n;
        }

        public static bool operator ==(EllipticCurveMontgomery E1, EllipticCurveMontgomery E2)
        {
            return E1.An == E2.An && E1.Ad == E2.Ad && E1.n == E2.n;
        }
        public static bool operator !=(EllipticCurveMontgomery E1, EllipticCurveMontgomery E2)
        {
            return !(E1 == E2);
        }
        public override bool Equals(Object obj)
        {
            if (obj == null)
                return false;
            EllipticCurveMontgomery E = obj as EllipticCurveMontgomery;
            return (new EllipticCurveMontgomery(this.An, this.Ad, this.n) == E);
        }
        public override int GetHashCode()
        {
            return base.GetHashCode();
        }
        public static EllipticCurveMontgomery GenerateCurveRandom(BigInteger n, out BigInteger x, out BigInteger y)
        {
            BigInteger a, b, d1, d2, inv;
            while (true)
            {
                a = Algorithms.RandomBigInteger(n) % n;
                if ((a + 2) % 4 != 0)
                    continue;
                x = Algorithms.RandomBigInteger(n) % n;
                y = Algorithms.RandomBigInteger(n) % n;
                d1 = (x * x * x + a * x * x + x) % n;
                d2 = (y * y) % n;
                inv = Algorithms.Invert(d2, n);
                if (inv < 0)
                    inv += n;
                if ((d2 * inv) % n != 1)
                    continue;
                b = (d1 * inv) % n;
                if ((b * (a * a - 4)) % n == 0)//A!={-2,2} & B!=0
                    continue;
                break;
            }
            return new EllipticCurveMontgomery(a, n);
        }
        public static EllipticCurveMontgomery GenerateCurve6(BigInteger n, out BigInteger x, out BigInteger z)// Z/6Z
        {
            BigInteger s, u, v, An, Ad, a, inv;
            var rand = new Random();
            while (true)
            {
                s = rand.Next() % (n - 6) + 6;
                u = (s * s - 5) % n;
                v = (4 * s) % n;
                BigInteger p = (u * u * u) % n;
                An = ((v - u) * (v - u) * (v - u) * (3 * u + v)) % n;
                Ad = (4 * p * v) % n;
                if (Ad < 0)
                    Ad += n;
                inv = Algorithms.Invert(Ad, n);
                if (inv * Ad % n != 1)
                    continue;
                a = (An * inv - 2) % n;
                x = p % n;
                z = (v * v * v) % n;

                inv = Algorithms.Invert(z, n);
                if (inv * z % n != 1)
                    continue;
                x = x * inv % n;
                break;
            }
            return new EllipticCurveMontgomery(a, n);
        }

    }
}
