using ECM2018.core.algorithms;
using ECM2018.core.weierstrass;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace ECM2018.core.edwards
{
    /// <summary>
    /// ax^2+y^2=1+d(x^2)(y^2)
    /// </summary>
    class EllipticCurveEdwards : IEllipticCurve
    {
        public BigInteger a { get; private set; }
        public BigInteger b { get; private set; }
        public BigInteger n { get; private set; }

        public EllipticCurveEdwards(BigInteger a, BigInteger d, BigInteger n)
        {
            this.a = a;
            this.b = d;
            this.n = n;

        }


        public static bool operator ==(EllipticCurveEdwards E1, EllipticCurveEdwards E2)
        {
            return E1.b == E2.b && E1.n == E2.n;
        }
        public static bool operator !=(EllipticCurveEdwards E1, EllipticCurveEdwards E2)
        {
            return !(E1 == E2);
        }
        public override bool Equals(Object obj)
        {
            if (obj == null)
                return false;
            EllipticCurveEdwards E = obj as EllipticCurveEdwards;
            return (new EllipticCurveEdwards(this.a, this.b, this.n) == E);
        }
        public override int GetHashCode()
        {
            return base.GetHashCode();
        }
        public static EllipticCurveEdwards GenerateCurveRandom(BigInteger n, out BigInteger x, out BigInteger y)
        {
            BigInteger a, d, d1, d2, inv;
            while (true)
            {
                x = Algorithms.RandomBigInteger(n) % n;
                y = Algorithms.RandomBigInteger(n) % n;
                if (x == 0 || y == 0)
                    continue;
                a = Algorithms.RandomBigInteger(n) % n;
                d1 = (a * x * x + y * y - 1) % n;
                d2 = (x * x * y * y) % n;
                inv = Algorithms.Invert(d2, n);
                if (inv < 0)
                    inv += n;
                if ((d2 * inv) % n != 1)
                    continue;
                d = (d1 * inv) % n;
                if (d == a)
                    continue;
                break;
            }
            return new EllipticCurveEdwards(a, d, n);
        }
        public static EllipticCurveEdwards GenerateCurve2x8(BigInteger n, out BigInteger x, out BigInteger y)// Z/2Z x Z/8Z
        {
            var E = new EllipticCurveWeierstrass(-8, -32, n);
            PointAffine P;
            BigInteger d;
            while (true)
            {
                P = new PointAffine(12, 40, E);
                BigInteger k = new Random().Next() % n;
                P = k * P;

                BigInteger temp = (P.x - 9) % n;
                if (temp < 0)
                    temp += n;
                BigInteger inv = Algorithms.Invert(temp, n);
                if (temp * inv % n != 1)
                    continue;
                temp = ((P.y + 25) * inv + 1) % n;
                if (temp < 0)
                    temp += n;
                BigInteger a = Algorithms.Invert(temp, n);
                if (temp * a % n != 1)
                    continue;

                temp = (8 * a * a - 1) % n;
                if (temp < 0)
                    temp += n;
                inv = Algorithms.Invert(temp, n);
                if (temp * inv % n != 1)
                    continue;
                BigInteger b = (2 * a * (4 * a + 1) * inv) % n;

                temp = BigInteger.ModPow(2 * b - 1, 4, n);
                if (temp < 0)
                    temp += n;
                inv = Algorithms.Invert(temp, n);
                if (temp * inv % n != 1)
                    continue;
                d = ((2 * (2 * b - 1) * (2 * b - 1) - 1) * inv) % n;

                temp = (6 * b - 5) % n;
                if (temp < 0)
                    temp += n;
                inv = Algorithms.Invert(temp, n);
                if (temp * inv % n != 1)
                    continue;
                x = ((2 * b - 1) * (4 * b - 3) * inv) % n;

                temp = ((P.y + 3 * P.x - 2) * (P.y + P.x + 16)) % n;
                if (temp < 0)
                    temp += n;
                inv = Algorithms.Invert(temp, n);
                if (temp * inv % n != 1)
                    continue;
                y = ((2 * b - 1) * (P.y * P.y + 50 * P.y - 2 * P.x * P.x * P.x + 27 * P.x * P.x - 104) * inv) % n;
                break;
            }

            return new EllipticCurveEdwards(1, d, n);
        }
        public static EllipticCurveEdwards GenerateCurve12(BigInteger n, out BigInteger x, out BigInteger y)// Z/12Z
        {
            var E = new EllipticCurveWeierstrass(-12, 0, n);
            PointAffine P;
            BigInteger d;
            while (true)
            {
                P = new PointAffine(6, -12, E);
                BigInteger k = new Random().Next() % n;
                P = k * P;
                if ((P.x == 0 && P.y == 0) || (P.x == -2 && BigInteger.Abs(P.y) == 4) || (P.x == 6 && BigInteger.Abs(P.y) == 12))
                    continue;

                BigInteger temp = (1024 * P.x * P.x * P.y * P.y) % n;
                BigInteger inv = Algorithms.Invert(temp, n);
                if (temp * inv % n != 1)
                    continue;
                temp = (-BigInteger.ModPow(P.x - 2, 3, n) * BigInteger.ModPow(P.x + 6, 3, n) * (P.x * P.x - 12 * P.x - 12)) % n;
                if (temp < 0)
                    temp += n;
                d = (temp * inv) % n;

                temp = ((P.x - 2) * (P.x + 6) * (P.x * P.x + 12 * P.x - 12)) % n;
                if (temp < 0)
                    temp += n;
                inv = Algorithms.Invert(temp, n);
                if (temp * inv % n != 1)
                    continue;
                x = (8 * P.y * (P.x * P.x + 12) * inv) % n;

                temp = ((P.x - 2) * (P.x + 6) * (P.x * P.x - 12)) % n;
                if (temp < 0)
                    temp += n;
                inv = Algorithms.Invert(temp, n);
                if (temp * inv % n != 1)
                    continue;
                y = (-4 * P.x * (P.x * P.x - 12 * P.x - 12) * inv) % n;

                break;
            }

            return new EllipticCurveEdwards(1, d, n);
        }
        public static EllipticCurveEdwards GetCurve(BigInteger n, out BigInteger x, out BigInteger y, string path)
        {
            string[] curves = System.IO.File.ReadAllLines(path);
            var rand = new Random();
            BigInteger temp, d, inv;
            while (true)
            {
                int num = rand.Next() % curves.Length;
                string[] param = curves[num].Split(' ');

                d = BigInteger.Parse(param[0]) % n;
                temp = BigInteger.Parse(param[1]) % n;
                inv = Algorithms.Invert(temp, n);
                if (inv < 0)
                    inv += n;
                if ((temp * inv) % n != 1)
                    continue;
                d = (d * inv) % n;

                x = BigInteger.Parse(param[2]) % n;
                temp = BigInteger.Parse(param[3]) % n;
                inv = Algorithms.Invert(temp, n);
                if (inv < 0)
                    inv += n;
                if ((temp * inv) % n != 1)
                {
                    continue;
                }
                x = (x * inv) % n;

                y = BigInteger.Parse(param[4]) % n;
                temp = BigInteger.Parse(param[5]) % n;
                inv = Algorithms.Invert(temp, n);
                if (inv < 0)
                    inv += n;
                if ((temp * inv) % n != 1)
                    continue;
                y = (y * inv) % n;

                break;
            }
            return new EllipticCurveEdwards(1, d, n);
        }

    }
}
