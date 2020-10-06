using ECM2018.core.algorithms;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace ECM2018.core.weierstrass
{
    class ECMWeierstrass : ECM
    {
        public static BigInteger GetFactorAffine(BigInteger n, int B1, int B2)
        {
            BigInteger x, y, g;
            var E = EllipticCurveWeierstrass.GenerateCurveRandom(n, out x, out y, out g, false);
            var P = new PointAffine(x, y, E);
            if (g > 1 && g < n)
            {
                return g;
            }
            //STEP 1
            var P1 = P as IPoint;
            var result = StageOneAffine(ref P1, n, B1);
            if (result > 1 || result == 0)
                return result;

            //STEP 2
            P = P1 as PointAffine;
            result = StageTwoAffine(ref P1, n, B1, B2);

            return result;
        }//first & second stage

        public static BigInteger GetFactorProjective(BigInteger n, int B1, int B2)
        {
            BigInteger x, y, g;
            var E = EllipticCurveWeierstrass.GenerateCurveRandom(n, out x, out y, out g, false);
            var P = new PointProjective(x, y, E);
            if (g > 1 && g < n)
            {
                return g;
            }

            //STEP 1
            var P1 = P as IPoint;
            var result = StageOneEdwardsProjective(ref P1, n, B1);
            if (result > 1 || result == 0)
                return result;

            //STEP 2
            result = StageTwoEdwardsProjective(ref P1, n, B1, B2);

            return result;
        }//projective

        public static BigInteger GetFactorJacobian(BigInteger n, int B1, int B2)
        {
            BigInteger x, y, g;
            var E = EllipticCurveWeierstrass.GenerateCurveRandom(n, out x, out y, out g, false);
            var P = new PointJacobian(x, y, E);
            if (g > 1 && g < n)
            {
                return g;
            }

            //STEP 1
            var P1 = P as IPoint;
            var result = StageOneProjective(ref P1, n, B1);
            if (result > 1 || result == 0)
                return result;

            //STEP 2
            result = StageTwoProjective(ref P1, n, B1, B2);

            return result;
        }//jacobian

    }
}
