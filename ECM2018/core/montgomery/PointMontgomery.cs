using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace ECM2018.core.montgomery
{
    /// <summary>
    /// B(Y^2)Z=X^3+A(X^2)Z+XZ^2
    /// </summary>
    class PointMontgomery : IPoint
    {
        public BigInteger x { get; private set; }
        public BigInteger y { get; private set; }
        public BigInteger z { get; private set; }

        public EllipticCurveMontgomery E { get; private set; }

        public PointMontgomery(BigInteger x, BigInteger y, BigInteger z, EllipticCurveMontgomery E)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            this.E = E;
        }

        public PointMontgomery(BigInteger x, BigInteger y, EllipticCurveMontgomery E)
        {
            this.x = x;
            this.y = y;
            this.z = 1;
            this.E = E;
        }

        public static bool operator ==(PointMontgomery P, PointMontgomery Q)
        {
            if (P.E == Q.E)
            {
                BigInteger v1 = (P.x * Q.z) % P.E.n;
                BigInteger v2 = (P.z * Q.x) % P.E.n;

                if (v1 < 0) v1 += P.E.n;
                if (v2 < 0) v2 += P.E.n;

                if (v1 == v2)
                        return true;
            }
            return false;
        }
        public static bool operator !=(PointMontgomery P, PointMontgomery Q)
        {
            return !(P == Q);
        }

        public static PointMontgomery operator *(BigInteger k, PointMontgomery P)
        {
            if (k == 0)
                return new PointMontgomery(0, 1, 0, P.E);
            else if (k == 1)
                return P;
            PointMontgomery U = P;
            PointMontgomery T = Double(P);
            if (k == 2)
                return T;
            var binary = new List<BigInteger>();
            while (k >= 2)
            {
                binary.Add(k % 2);
                k /= 2;
            }
            PointMontgomery R = P;
            for (int i = binary.Count - 1; i >= 0; i--)
            {
                if (binary[i] == 1)
                {
                    U = Addition(T, U, P);
                    T = Double(T);
                }
                else
                {
                    T = Addition(T, U, P);
                    U = Double(U);
                }
            }
            return U;
        }

        public static PointMontgomery Double(PointMontgomery P)//2M + 2S + 1*a24 + 4add
        {
            BigInteger a = P.x + P.z;
            BigInteger b = P.x - P.z;
            BigInteger t1 = a * a;
            BigInteger t2 = b * b;
            BigInteger t = t1 - t2;
            BigInteger x = (t1 * t2) % P.E.n;
            BigInteger z = (t * (t2 + t * P.E.b)) % P.E.n;
            return new PointMontgomery(x, 0, z, P.E);
        }

        public static PointMontgomery Addition(PointMontgomery P, PointMontgomery Q, PointMontgomery P0)//4M + 2S + 6add
        {
            BigInteger t1 = (P.x - P.z) * (Q.x + Q.z);
            BigInteger t2 = (P.x + P.z) * (Q.x - Q.z);
            BigInteger a = t1 + t2;
            BigInteger b = t1 - t2;
            BigInteger x = (P0.z * a * a) % P.E.n;
            BigInteger z = (P0.x * b * b) % P.E.n;
            return new PointMontgomery(x, 0, z, P.E);
        }

        public bool IsOnCurve()
        {
            return ((x * x * x + E.a * x * x * z + x * z * z) % E.n == 0);
        }


        public override bool Equals(Object obj)
        {
            throw new NotImplementedException();
        }
        public override int GetHashCode()
        {
            return this.GetHashCode();
        }
        public bool IsInfinite()
        {
            return x == 0 && z == 0;
        }
        public IPoint Multiply(BigInteger n)
        {
            return n * (this as PointMontgomery);
        }
        public IPoint Add(IPoint P)
        {
            throw new NotImplementedException();
        }
        public IPoint Double()
        {
            return this.Double();
        }
        public IPoint Clone()
        {
            return new PointMontgomery(this.x, this.y, this.z, this.E);
        }
        public void ToPointAffine()
        {
        }

    }
}

