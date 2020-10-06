using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace ECM2018.core.montgomery
{
    /// <summary>
    /// x^3+(An/Ad-2)(x^2)z+xz^2(mod n)=0
    /// </summary>
    class PointMontgomeryBrent : IPoint
    {
        public BigInteger x { get; private set; }
        public BigInteger y { get; private set; }
        public BigInteger z { get; private set; }

        public IEllipticCurve E { get; private set; }

        public PointMontgomeryBrent(BigInteger x, BigInteger y, BigInteger z, IEllipticCurve E)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            this.E = E;
        }

        public static bool operator ==(PointMontgomeryBrent P, PointMontgomeryBrent Q)//todo: доделать
        {
            if (P.E == Q.E && P.x == Q.x && P.y == Q.y)
            {
                return true;
            }
            return false;
        }
        public static bool operator !=(PointMontgomeryBrent P, PointMontgomeryBrent Q)
        {
            return !(P == Q);
        }

        public static PointMontgomeryBrent operator *(BigInteger k, PointMontgomeryBrent P)
        {
            if (k == 0)
                return new PointMontgomeryBrent(0, 1, 0, P.E);
            else if (k == 1)
                return P;
            PointMontgomeryBrent U = P;
            PointMontgomeryBrent T = Double(P);
            if (k == 2)
                return T;
            List<BigInteger> binary = new List<BigInteger>();
            while (k >= 2)
            {
                binary.Add(k % 2);
                k /= 2;
            }
            PointMontgomeryBrent R = P;
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

        public static PointMontgomeryBrent Double(PointMontgomeryBrent P)
        {
            EllipticCurveMontgomery E = P.E as EllipticCurveMontgomery;
            BigInteger t1 = (P.x + P.z) * (P.x + P.z);
            BigInteger t2 = (P.x - P.z) * (P.x - P.z);
            BigInteger t = t1 - t2;
            BigInteger x = (t1 * t2 * 4 * E.Ad) % P.E.n;
            BigInteger z = ((4 * E.Ad * t2 + t * E.An) * t) % P.E.n;
            return new PointMontgomeryBrent(x, 0, z, P.E);
        }

        public static PointMontgomeryBrent Addition(PointMontgomeryBrent P, PointMontgomeryBrent Q, PointMontgomeryBrent P0)
        {
            BigInteger t1 = (P.x - P.z) * (Q.x + Q.z);
            BigInteger t2 = (P.x + P.z) * (Q.x - Q.z);
            BigInteger x = (P0.z * (t1 + t2) * (t1 + t2)) % P.E.n;
            BigInteger z = (P0.x * (t1 - t2) * (t1 - t2)) % P.E.n;
            return new PointMontgomeryBrent(x, 0, z, P.E);
        }

        public bool IsOnCurve()
        {
            EllipticCurveMontgomery E = this.E as EllipticCurveMontgomery;
            return ((x * x * x + (E.An / E.Ad - 2) * x * x * z + x * z * z) % E.n == 0);
        }


        public override bool Equals(Object obj)
        {
            if (obj == null)
                return false;
            PointMontgomeryBrent P = obj as PointMontgomeryBrent;
            return (new PointMontgomeryBrent(this.x, this.y, this.z, this.E) == P);
        }
        public override int GetHashCode()
        {
            return base.GetHashCode();
        }
        public bool IsInfinite()
        {
            return x == 0 && z == 0;
        }
        public IPoint Multiply(BigInteger n)
        {
            return n * (this as PointMontgomeryBrent);
        }
        public IPoint Add(IPoint P)
        {
            return this;
        }
        public IPoint Double()
        {
            return this.Double();
        }
    }
}
