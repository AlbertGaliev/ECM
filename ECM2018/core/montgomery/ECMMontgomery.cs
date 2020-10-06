using ECM2018.core.algorithms;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace ECM2018.core.montgomery
{
    class ECMMontgomery : ECM
    {
        public static BigInteger GetFactorProjective(BigInteger n, int B1, int B2, TorsionType type)
        {
            BigInteger x;
            BigInteger y;
            var E = GenerateCurve(n, out x, out y, type);
            var P = new PointMontgomery(x, y, E);
            //STEP 1
            var P1 = P as IPoint;
            var result = StageOneProjective(ref P1, n, B1);
            if (result > 1 || result == 0)
                return result;

            //STEP 2
            result = StageTwoProjective2(ref P1, n, B1, B2);

            return result;

        }//Montgomery

        private static EllipticCurveMontgomery GenerateCurve(BigInteger n, out BigInteger x, out BigInteger y, TorsionType type)
        {
            switch (type)
            {
                case TorsionType.RandomCurve:
                    return EllipticCurveMontgomery.GenerateCurveRandom(n, out x, out y);
                case TorsionType.Curve6:
                    return EllipticCurveMontgomery.GenerateCurve6(n, out x, out y);
                default:
                    return EllipticCurveMontgomery.GenerateCurveRandom(n, out x, out y);
            }
        }

    }
}
