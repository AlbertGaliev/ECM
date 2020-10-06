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
    class ECMEdwards : ECM
    {
        public static BigInteger GetFactorProjective(BigInteger n, int B1, int B2, TorsionType type)
        {
            BigInteger x, y;
            var E = GenerateCurve(n, out x, out y, type);
            var P = new PointEdwardsProjective(x, y, E);

            //STEP 1
            var P1 = P as IPoint;
            var result = StageOneEdwardsProjective(ref P1, n, B1);
            if (result > 1 || result == 0)
                return result;

            //STEP 2
            result = StageTwoEdwardsProjective(ref P1, n, B1, B2);

            return result;
        }//projective

        public static BigInteger GetFactorInverted(BigInteger n, int B1, int B2, TorsionType type)
        {
            BigInteger x, y;
            var E = GenerateCurve(n, out x, out y, type);
            var P = new PointEdwardsInverted(x, y, E);

            //STEP 1
            var P1 = P as IPoint;
            var result = StageOneEdwardsProjective(ref P1, n, B1);
            if (result > 1 || result == 0)
                return result;

            //STEP 2
            result = StageTwoEdwardsProjective(ref P1, n, B1, B2);

            return result;
        }//inverted

        public static BigInteger GetFactorExtended(BigInteger n, int B1, int B2, TorsionType type)
        {
            BigInteger x, y;
            var E = GenerateCurve(n, out x, out y, type);
            var P = new PointEdwardsExtended(x, y, E);

            //STEP 1
            var P1 = P as IPoint;
            var result = StageOneEdwardsProjective(ref P1, n, B1);
            if (result > 1 || result == 0)
                return result;

            //STEP 2
            result = StageTwoEdwardsProjective(ref P1, n, B1, B2);

            return result;
        }//extended

        private static EllipticCurveEdwards GenerateCurve(BigInteger n, out BigInteger x, out BigInteger y, TorsionType type)
        {
            switch (type)
            {
                case TorsionType.RandomCurve:
                    return EllipticCurveEdwards.GenerateCurveRandom(n, out x, out y);
                case TorsionType.Curve12:
                    return EllipticCurveEdwards.GenerateCurve12(n, out x, out y);
                case TorsionType.Curve2x8:
                    return EllipticCurveEdwards.GenerateCurve2x8(n, out x, out y);
                case TorsionType.Curve12FromFile:
                    return EllipticCurveEdwards.GetCurve(n, out x, out y, "curve12.txt");
                case TorsionType.Curve2x8FromFile:
                    return EllipticCurveEdwards.GetCurve(n, out x, out y, "curve2x8.txt");
                default:
                    return EllipticCurveEdwards.GenerateCurveRandom(n, out x, out y);
            }
        }
    }
}
