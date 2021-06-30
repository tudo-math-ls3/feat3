// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_GAUSS_LEGENDRE_DRIVER_HPP
#define KERNEL_CUBATURE_SCALAR_GAUSS_LEGENDRE_DRIVER_HPP 1

// includes, FEAT
#include <kernel/cubature/scalar/driver_base.hpp>
#include <kernel/util/math.hpp>

namespace FEAT
{
  namespace Cubature
  {
    namespace Scalar
    {
      /**
       * \brief Gauss-Legendre Rule driver class template
       *
       * This driver implements the Gauss-Legendre rule.
       * \see http://en.wikipedia.org/wiki/Gaussian_quadrature
       * \see http://mathworld.wolfram.com/GaussianQuadrature.html
       *
       * \author Peter Zajac
       */
      class GaussLegendreDriver :
        public DriverBase
      {
      public:
        /// this rule is variadic
        static constexpr bool variadic = true;
        /// this rule has at least 1 point
        static constexpr int min_points = 1;
        /// this rule has at most 5 points
        static constexpr int max_points = 20;

        /// Returns the name of the cubature rule.
        static String name()
        {
          return "gauss-legendre";
        }

        /**
         * \brief Fills the cubature rule structure.
         *
         * \param[in,out] rule
         * The cubature rule to be filled.
         *
         * \param[in] num_points
         * The number of quadrature points.
         */
        template<
          typename Weight_,
          typename Coord_>
        static void fill(Rule<Weight_, Coord_>& rule, int num_points)
        {
          // auxiliary variables
          //Coord_ dc;
          //Weight_ dw;

          // how many points do we have?
          switch(num_points)
          {
          /*case 1:
            rule.get_coord(0) = Coord_(0);

            rule.get_weight(0) = Weight_(2);
            break;

          case 2:
            rule.get_coord(0) = -Math::sqrt(Coord_(1) / Coord_(3));
            rule.get_coord(1) = +Math::sqrt(Coord_(1) / Coord_(3));

            rule.get_weight(0) = Weight_(1);
            rule.get_weight(1) = Weight_(1);
            break;

          case 3:
            rule.get_coord(0) = -Math::sqrt(Coord_(3) / Coord_(5));
            rule.get_coord(1) = Coord_(0);
            rule.get_coord(2) = +Math::sqrt(Coord_(3) / Coord_(5));

            rule.get_weight(0) = Weight_(5) / Weight_(9);
            rule.get_weight(1) = Weight_(8) / Weight_(9);
            rule.get_weight(2) = Weight_(5) / Weight_(9);
            break;

          case 4:
            dc = Math::sqrt(Coord_(24) / Coord_(5));
            rule.get_coord(0) = -Math::sqrt((Coord_(3) + dc) / Coord_(7));
            rule.get_coord(1) = -Math::sqrt((Coord_(3) - dc) / Coord_(7));
            rule.get_coord(2) = +Math::sqrt((Coord_(3) - dc) / Coord_(7));
            rule.get_coord(3) = +Math::sqrt((Coord_(3) + dc) / Coord_(7));

            dw = Math::sqrt(Weight_(30));
            rule.get_weight(0) = (Weight_(18) - dw) / Weight_(36);
            rule.get_weight(1) = (Weight_(18) + dw) / Weight_(36);
            rule.get_weight(2) = (Weight_(18) + dw) / Weight_(36);
            rule.get_weight(3) = (Weight_(18) - dw) / Weight_(36);
            break;

          case 5:
            dc = Coord_(2) * Math::sqrt(Coord_(10) / Coord_(7));
            rule.get_coord(0) = -Math::sqrt(Coord_(5) + dc) / Coord_(3);
            rule.get_coord(1) = -Math::sqrt(Coord_(5) - dc) / Coord_(3);
            rule.get_coord(2) = Coord_(0);
            rule.get_coord(3) = +Math::sqrt(Coord_(5) - dc) / Coord_(3);
            rule.get_coord(4) = +Math::sqrt(Coord_(5) + dc) / Coord_(3);

            dw = Weight_(13) * Math::sqrt(Weight_(70));
            rule.get_weight(0) = (Weight_(322) - dw) / Weight_(900);
            rule.get_weight(1) = (Weight_(322) + dw) / Weight_(900);
            rule.get_weight(2) =  Weight_(128)       / Weight_(225);
            rule.get_weight(3) = (Weight_(322) + dw) / Weight_(900);
            rule.get_weight(4) = (Weight_(322) - dw) / Weight_(900);
            break;*/
          case 1:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(+0.0));
            rule.get_weight( 0) = Weight_(FEAT_F128C(2.0));
            break;

          case 2:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(-0.57735026918962576450914878050195745564760175127));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(+0.57735026918962576450914878050195745564760175127));
            rule.get_weight( 0) = Weight_(FEAT_F128C(1.0));
            rule.get_weight( 1) = Weight_(FEAT_F128C(1.0));
            break;

          case 3:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(+0.0));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(-0.77459666924148337703585307995647992216658434105));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(+0.77459666924148337703585307995647992216658434105));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.888888888888888888888888888888888888888888888888));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.555555555555555555555555555555555555555555555555));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.555555555555555555555555555555555555555555555555));
            break;

          case 4:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(-0.33998104358485626480266575910324468720057586977));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(+0.33998104358485626480266575910324468720057586977));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(-0.86113631159405257522394648889280950509572537962));
            rule.get_coord ( 3) = Coord_(FEAT_F128C(+0.86113631159405257522394648889280950509572537962));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.652145154862546142626936050778000592764651304166));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.652145154862546142626936050778000592764651304166));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.347854845137453857373063949221999407235348695833));
            rule.get_weight( 3) = Weight_(FEAT_F128C(0.347854845137453857373063949221999407235348695833));
            break;

          case 5:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(+0.0));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(-0.53846931010568309103631442070020880496728660690));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(+0.53846931010568309103631442070020880496728660690));
            rule.get_coord ( 3) = Coord_(FEAT_F128C(-0.90617984593866399279762687829939296512565191076));
            rule.get_coord ( 4) = Coord_(FEAT_F128C(+0.90617984593866399279762687829939296512565191076));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.568888888888888888888888888888888888888888888888));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.478628670499366468041291514835638192912295553343));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.478628670499366468041291514835638192912295553343));
            rule.get_weight( 3) = Weight_(FEAT_F128C(0.236926885056189087514264040719917362643260002212));
            rule.get_weight( 4) = Weight_(FEAT_F128C(0.236926885056189087514264040719917362643260002212));
            break;

          case 6:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(+0.66120938646626451366139959501990534700644856439));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(-0.66120938646626451366139959501990534700644856439));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(-0.23861918608319690863050172168071193541861063014));
            rule.get_coord ( 3) = Coord_(FEAT_F128C(+0.23861918608319690863050172168071193541861063014));
            rule.get_coord ( 4) = Coord_(FEAT_F128C(-0.93246951420315202781230155449399460913476573771));
            rule.get_coord ( 5) = Coord_(FEAT_F128C(+0.93246951420315202781230155449399460913476573771));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.360761573048138607569833513837716111661521892746));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.360761573048138607569833513837716111661521892746));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.467913934572691047389870343989550994811655605769));
            rule.get_weight( 3) = Weight_(FEAT_F128C(0.467913934572691047389870343989550994811655605769));
            rule.get_weight( 4) = Weight_(FEAT_F128C(0.171324492379170345040296142172732893526822501484));
            rule.get_weight( 5) = Weight_(FEAT_F128C(0.171324492379170345040296142172732893526822501484));
            break;

          case 7:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(+0.0));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(+0.40584515137739716690660641207696146334738201409));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(-0.40584515137739716690660641207696146334738201409));
            rule.get_coord ( 3) = Coord_(FEAT_F128C(-0.74153118559939443986386477328078840707414764714));
            rule.get_coord ( 4) = Coord_(FEAT_F128C(+0.74153118559939443986386477328078840707414764714));
            rule.get_coord ( 5) = Coord_(FEAT_F128C(-0.94910791234275852452618968404785126240077093767));
            rule.get_coord ( 6) = Coord_(FEAT_F128C(+0.94910791234275852452618968404785126240077093767));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.417959183673469387755102040816326530612244897959));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.381830050505118944950369775488975133878365083533));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.381830050505118944950369775488975133878365083533));
            rule.get_weight( 3) = Weight_(FEAT_F128C(0.279705391489276667901467771423779582486925065226));
            rule.get_weight( 4) = Weight_(FEAT_F128C(0.279705391489276667901467771423779582486925065226));
            rule.get_weight( 5) = Weight_(FEAT_F128C(0.129484966168869693270611432679082018328587402259));
            rule.get_weight( 6) = Weight_(FEAT_F128C(0.129484966168869693270611432679082018328587402259));
            break;

          case 8:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(-0.18343464249564980493947614236018398066675781291));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(+0.18343464249564980493947614236018398066675781291));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(-0.52553240991632898581773904918924634904196424312));
            rule.get_coord ( 3) = Coord_(FEAT_F128C(+0.52553240991632898581773904918924634904196424312));
            rule.get_coord ( 4) = Coord_(FEAT_F128C(-0.79666647741362673959155393647583043683717173161));
            rule.get_coord ( 5) = Coord_(FEAT_F128C(+0.79666647741362673959155393647583043683717173161));
            rule.get_coord ( 6) = Coord_(FEAT_F128C(-0.96028985649753623168356086856947299042823523430));
            rule.get_coord ( 7) = Coord_(FEAT_F128C(+0.96028985649753623168356086856947299042823523430));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.362683783378361982965150449277195612194146039894));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.362683783378361982965150449277195612194146039894));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.313706645877887287337962201986601313260328999002));
            rule.get_weight( 3) = Weight_(FEAT_F128C(0.313706645877887287337962201986601313260328999002));
            rule.get_weight( 4) = Weight_(FEAT_F128C(0.222381034453374470544355994426240884430130870051));
            rule.get_weight( 5) = Weight_(FEAT_F128C(0.222381034453374470544355994426240884430130870051));
            rule.get_weight( 6) = Weight_(FEAT_F128C(0.101228536290376259152531354309962190115394091051));
            rule.get_weight( 7) = Weight_(FEAT_F128C(0.101228536290376259152531354309962190115394091051));
            break;

          case 9:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(+0.0));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(-0.83603110732663579429942978806973487654410671812));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(+0.83603110732663579429942978806973487654410671812));
            rule.get_coord ( 3) = Coord_(FEAT_F128C(-0.96816023950762608983557620290367287004940480049));
            rule.get_coord ( 4) = Coord_(FEAT_F128C(+0.96816023950762608983557620290367287004940480049));
            rule.get_coord ( 5) = Coord_(FEAT_F128C(-0.32425342340380892903853801464333660857195626073));
            rule.get_coord ( 6) = Coord_(FEAT_F128C(+0.32425342340380892903853801464333660857195626073));
            rule.get_coord ( 7) = Coord_(FEAT_F128C(-0.61337143270059039730870203934147418478572060494));
            rule.get_coord ( 8) = Coord_(FEAT_F128C(+0.61337143270059039730870203934147418478572060494));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.330239355001259763164525069286974048878810783572));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.180648160694857404058472031242912809514337821732));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.180648160694857404058472031242912809514337821732));
            rule.get_weight( 3) = Weight_(FEAT_F128C(0.081274388361574411971892158110523650675661720782));
            rule.get_weight( 4) = Weight_(FEAT_F128C(0.081274388361574411971892158110523650675661720782));
            rule.get_weight( 5) = Weight_(FEAT_F128C(0.312347077040002840068630406584443665598754861261));
            rule.get_weight( 6) = Weight_(FEAT_F128C(0.312347077040002840068630406584443665598754861261));
            rule.get_weight( 7) = Weight_(FEAT_F128C(0.260610696402935462318742869418632849771840204437));
            rule.get_weight( 8) = Weight_(FEAT_F128C(0.260610696402935462318742869418632849771840204437));
            break;

          case 10:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(-0.14887433898163121088482600112971998461756485942));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(+0.14887433898163121088482600112971998461756485942));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(-0.43339539412924719079926594316578416220007183765));
            rule.get_coord ( 3) = Coord_(FEAT_F128C(+0.43339539412924719079926594316578416220007183765));
            rule.get_coord ( 4) = Coord_(FEAT_F128C(-0.67940956829902440623432736511487357576929471183));
            rule.get_coord ( 5) = Coord_(FEAT_F128C(+0.67940956829902440623432736511487357576929471183));
            rule.get_coord ( 6) = Coord_(FEAT_F128C(-0.86506336668898451073209668842349304852754301496));
            rule.get_coord ( 7) = Coord_(FEAT_F128C(+0.86506336668898451073209668842349304852754301496));
            rule.get_coord ( 8) = Coord_(FEAT_F128C(-0.97390652851717172007796401208445205342826994669));
            rule.get_coord ( 9) = Coord_(FEAT_F128C(+0.97390652851717172007796401208445205342826994669));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.295524224714752870173892994651338329421046717026));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.295524224714752870173892994651338329421046717026));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.269266719309996355091226921569469352859759938460));
            rule.get_weight( 3) = Weight_(FEAT_F128C(0.269266719309996355091226921569469352859759938460));
            rule.get_weight( 4) = Weight_(FEAT_F128C(0.219086362515982043995534934228163192458771870522));
            rule.get_weight( 5) = Weight_(FEAT_F128C(0.219086362515982043995534934228163192458771870522));
            rule.get_weight( 6) = Weight_(FEAT_F128C(0.149451349150580593145776339657697332402556639669));
            rule.get_weight( 7) = Weight_(FEAT_F128C(0.149451349150580593145776339657697332402556639669));
            rule.get_weight( 8) = Weight_(FEAT_F128C(0.066671344308688137593568809893331792857864834320));
            rule.get_weight( 9) = Weight_(FEAT_F128C(0.066671344308688137593568809893331792857864834320));
            break;

          case 11:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(+0.0));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(-0.26954315595234497233153198540086152467962186243));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(+0.26954315595234497233153198540086152467962186243));
            rule.get_coord ( 3) = Coord_(FEAT_F128C(-0.51909612920681181592572566945860955448022711511));
            rule.get_coord ( 4) = Coord_(FEAT_F128C(+0.51909612920681181592572566945860955448022711511));
            rule.get_coord ( 5) = Coord_(FEAT_F128C(-0.73015200557404932409341625203115345804964306202));
            rule.get_coord ( 6) = Coord_(FEAT_F128C(+0.73015200557404932409341625203115345804964306202));
            rule.get_coord ( 7) = Coord_(FEAT_F128C(-0.88706259976809529907515776930392726663167575122));
            rule.get_coord ( 8) = Coord_(FEAT_F128C(+0.88706259976809529907515776930392726663167575122));
            rule.get_coord ( 9) = Coord_(FEAT_F128C(-0.97822865814605699280393800112285739077142240891));
            rule.get_coord (10) = Coord_(FEAT_F128C(+0.97822865814605699280393800112285739077142240891));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.272925086777900630714483528336342189156041969894));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.262804544510246662180688869890509195372764677603));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.262804544510246662180688869890509195372764677603));
            rule.get_weight( 3) = Weight_(FEAT_F128C(0.233193764591990479918523704843175139431798172316));
            rule.get_weight( 4) = Weight_(FEAT_F128C(0.233193764591990479918523704843175139431798172316));
            rule.get_weight( 5) = Weight_(FEAT_F128C(0.186290210927734251426097641431655891691284748040));
            rule.get_weight( 6) = Weight_(FEAT_F128C(0.186290210927734251426097641431655891691284748040));
            rule.get_weight( 7) = Weight_(FEAT_F128C(0.125580369464904624634694299223940100197615791395));
            rule.get_weight( 8) = Weight_(FEAT_F128C(0.125580369464904624634694299223940100197615791395));
            rule.get_weight( 9) = Weight_(FEAT_F128C(0.055668567116173666482753720442548578728515625696));
            rule.get_weight(10) = Weight_(FEAT_F128C(0.055668567116173666482753720442548578728515625696));
            break;

          case 12:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(-0.12523340851146891547244136946385312998339691630));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(+0.12523340851146891547244136946385312998339691630));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(-0.36783149899818019375269153664371756125636014133));
            rule.get_coord ( 3) = Coord_(FEAT_F128C(+0.36783149899818019375269153664371756125636014133));
            rule.get_coord ( 4) = Coord_(FEAT_F128C(-0.58731795428661744729670241894053428036909851404));
            rule.get_coord ( 5) = Coord_(FEAT_F128C(+0.58731795428661744729670241894053428036909851404));
            rule.get_coord ( 6) = Coord_(FEAT_F128C(-0.76990267419430468703689383321281807598492575001));
            rule.get_coord ( 7) = Coord_(FEAT_F128C(+0.76990267419430468703689383321281807598492575001));
            rule.get_coord ( 8) = Coord_(FEAT_F128C(-0.90411725637047485667846586611909619253759670921));
            rule.get_coord ( 9) = Coord_(FEAT_F128C(+0.90411725637047485667846586611909619253759670921));
            rule.get_coord (10) = Coord_(FEAT_F128C(-0.98156063424671925069054909014928082296015519981));
            rule.get_coord (11) = Coord_(FEAT_F128C(+0.98156063424671925069054909014928082296015519981));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.249147045813402785000562436042951210830460902569));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.249147045813402785000562436042951210830460902569));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.233492536538354808760849898924878056259409972199));
            rule.get_weight( 3) = Weight_(FEAT_F128C(0.233492536538354808760849898924878056259409972199));
            rule.get_weight( 4) = Weight_(FEAT_F128C(0.203167426723065921749064455809798376506518147274));
            rule.get_weight( 5) = Weight_(FEAT_F128C(0.203167426723065921749064455809798376506518147274));
            rule.get_weight( 6) = Weight_(FEAT_F128C(0.160078328543346226334652529543359071872011730490));
            rule.get_weight( 7) = Weight_(FEAT_F128C(0.160078328543346226334652529543359071872011730490));
            rule.get_weight( 8) = Weight_(FEAT_F128C(0.106939325995318430960254718193996224214570173470));
            rule.get_weight( 9) = Weight_(FEAT_F128C(0.106939325995318430960254718193996224214570173470));
            rule.get_weight(10) = Weight_(FEAT_F128C(0.047175336386511827194615961485017060317029073994));
            rule.get_weight(11) = Weight_(FEAT_F128C(0.047175336386511827194615961485017060317029073994));
            break;

          case 13:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(+0.0));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(-0.23045831595513479406552812109798883521154237588));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(+0.23045831595513479406552812109798883521154237588));
            rule.get_coord ( 3) = Coord_(FEAT_F128C(-0.44849275103644685287791285212763986780192166744));
            rule.get_coord ( 4) = Coord_(FEAT_F128C(+0.44849275103644685287791285212763986780192166744));
            rule.get_coord ( 5) = Coord_(FEAT_F128C(-0.64234933944034022064398460699551565007169739826));
            rule.get_coord ( 6) = Coord_(FEAT_F128C(+0.64234933944034022064398460699551565007169739826));
            rule.get_coord ( 7) = Coord_(FEAT_F128C(-0.80157809073330991279420648958285989030561572479));
            rule.get_coord ( 8) = Coord_(FEAT_F128C(+0.80157809073330991279420648958285989030561572479));
            rule.get_coord ( 9) = Coord_(FEAT_F128C(-0.91759839922297796520654783650071951239047479011));
            rule.get_coord (10) = Coord_(FEAT_F128C(+0.91759839922297796520654783650071951239047479011));
            rule.get_coord (11) = Coord_(FEAT_F128C(-0.98418305471858814947282944880710961106499056192));
            rule.get_coord (12) = Coord_(FEAT_F128C(+0.98418305471858814947282944880710961106499056192));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.232551553230873910194589515268835948156627477306));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.226283180262897238412090186039776618434757737615));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.226283180262897238412090186039776618434757737615));
            rule.get_weight( 3) = Weight_(FEAT_F128C(0.207816047536888502312523219306052763386582609199));
            rule.get_weight( 4) = Weight_(FEAT_F128C(0.207816047536888502312523219306052763386582609199));
            rule.get_weight( 5) = Weight_(FEAT_F128C(0.178145980761945738280046691996097995512812650661));
            rule.get_weight( 6) = Weight_(FEAT_F128C(0.178145980761945738280046691996097995512812650661));
            rule.get_weight( 7) = Weight_(FEAT_F128C(0.138873510219787238463601776868871467621862718263));
            rule.get_weight( 8) = Weight_(FEAT_F128C(0.138873510219787238463601776868871467621862718263));
            rule.get_weight( 9) = Weight_(FEAT_F128C(0.092121499837728447914421775953797120923683999862));
            rule.get_weight(10) = Weight_(FEAT_F128C(0.092121499837728447914421775953797120923683999862));
            rule.get_weight(11) = Weight_(FEAT_F128C(0.040484004765315879520021592200986060041986545744));
            rule.get_weight(12) = Weight_(FEAT_F128C(0.040484004765315879520021592200986060041986545744));
            break;

          case 14:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(-0.10805494870734366206624465021983474761195160547));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(+0.10805494870734366206624465021983474761195160547));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(-0.31911236892788976043567182416847546683426120353));
            rule.get_coord ( 3) = Coord_(FEAT_F128C(+0.31911236892788976043567182416847546683426120353));
            rule.get_coord ( 4) = Coord_(FEAT_F128C(-0.51524863635815409196529071855118866230888528256));
            rule.get_coord ( 5) = Coord_(FEAT_F128C(+0.51524863635815409196529071855118866230888528256));
            rule.get_coord ( 6) = Coord_(FEAT_F128C(-0.68729290481168547014801980301933413753840121274));
            rule.get_coord ( 7) = Coord_(FEAT_F128C(+0.68729290481168547014801980301933413753840121274));
            rule.get_coord ( 8) = Coord_(FEAT_F128C(-0.82720131506976499318979474265039496103970110147));
            rule.get_coord ( 9) = Coord_(FEAT_F128C(+0.82720131506976499318979474265039496103970110147));
            rule.get_coord (10) = Coord_(FEAT_F128C(-0.92843488366357351733639113937787426447703921040));
            rule.get_coord (11) = Coord_(FEAT_F128C(+0.92843488366357351733639113937787426447703921040));
            rule.get_coord (12) = Coord_(FEAT_F128C(-0.98628380869681233884159726670405280167609140723));
            rule.get_coord (13) = Coord_(FEAT_F128C(+0.98628380869681233884159726670405280167609140723));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.215263853463157790195876443316260035274997558054));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.215263853463157790195876443316260035274997558054));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.205198463721295603965924065661218055710339061309));
            rule.get_weight( 3) = Weight_(FEAT_F128C(0.205198463721295603965924065661218055710339061309));
            rule.get_weight( 4) = Weight_(FEAT_F128C(0.185538397477937813741716590125157036248922602937));
            rule.get_weight( 5) = Weight_(FEAT_F128C(0.185538397477937813741716590125157036248922602937));
            rule.get_weight( 6) = Weight_(FEAT_F128C(0.157203167158193534569601938623842156605668037337));
            rule.get_weight( 7) = Weight_(FEAT_F128C(0.157203167158193534569601938623842156605668037337));
            rule.get_weight( 8) = Weight_(FEAT_F128C(0.121518570687903184689414809072476625956669345690));
            rule.get_weight( 9) = Weight_(FEAT_F128C(0.121518570687903184689414809072476625956669345690));
            rule.get_weight(10) = Weight_(FEAT_F128C(0.080158087159760209805633277062854309583697785394));
            rule.get_weight(11) = Weight_(FEAT_F128C(0.080158087159760209805633277062854309583697785394));
            rule.get_weight(12) = Weight_(FEAT_F128C(0.035119460331751863031832876138191780619705609277));
            rule.get_weight(13) = Weight_(FEAT_F128C(0.035119460331751863031832876138191780619705609277));
            break;

          case 15:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(+0.0));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(-0.20119409399743452230062830339459620781283645446));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(+0.20119409399743452230062830339459620781283645446));
            rule.get_coord ( 3) = Coord_(FEAT_F128C(-0.39415134707756336989720737098104546836275277615));
            rule.get_coord ( 4) = Coord_(FEAT_F128C(+0.39415134707756336989720737098104546836275277615));
            rule.get_coord ( 5) = Coord_(FEAT_F128C(-0.57097217260853884753722673725391064123838639628));
            rule.get_coord ( 6) = Coord_(FEAT_F128C(+0.57097217260853884753722673725391064123838639628));
            rule.get_coord ( 7) = Coord_(FEAT_F128C(-0.72441773136017004741618605461393800963089929458));
            rule.get_coord ( 8) = Coord_(FEAT_F128C(+0.72441773136017004741618605461393800963089929458));
            rule.get_coord ( 9) = Coord_(FEAT_F128C(-0.84820658341042721620064832077421685136625617473));
            rule.get_coord (10) = Coord_(FEAT_F128C(+0.84820658341042721620064832077421685136625617473));
            rule.get_coord (11) = Coord_(FEAT_F128C(-0.93727339240070590430775894771020947124399627351));
            rule.get_coord (12) = Coord_(FEAT_F128C(+0.93727339240070590430775894771020947124399627351));
            rule.get_coord (13) = Coord_(FEAT_F128C(-0.98799251802048542848956571858661258114697281712));
            rule.get_coord (14) = Coord_(FEAT_F128C(+0.98799251802048542848956571858661258114697281712));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.202578241925561272880620199967519314838662158009));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.198431485327111576456118326443839324818692559957));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.198431485327111576456118326443839324818692559957));
            rule.get_weight( 3) = Weight_(FEAT_F128C(0.186161000015562211026800561866422824506226012277));
            rule.get_weight( 4) = Weight_(FEAT_F128C(0.186161000015562211026800561866422824506226012277));
            rule.get_weight( 5) = Weight_(FEAT_F128C(0.166269205816993933553200860481208811130900180098));
            rule.get_weight( 6) = Weight_(FEAT_F128C(0.166269205816993933553200860481208811130900180098));
            rule.get_weight( 7) = Weight_(FEAT_F128C(0.139570677926154314447804794511028322520850275315));
            rule.get_weight( 8) = Weight_(FEAT_F128C(0.139570677926154314447804794511028322520850275315));
            rule.get_weight( 9) = Weight_(FEAT_F128C(0.107159220467171935011869546685869303415543715758));
            rule.get_weight(10) = Weight_(FEAT_F128C(0.107159220467171935011869546685869303415543715758));
            rule.get_weight(11) = Weight_(FEAT_F128C(0.070366047488108124709267416450667338466708032754));
            rule.get_weight(12) = Weight_(FEAT_F128C(0.070366047488108124709267416450667338466708032754));
            rule.get_weight(13) = Weight_(FEAT_F128C(0.030753241996117268354628393577204417721748144833));
            rule.get_weight(14) = Weight_(FEAT_F128C(0.030753241996117268354628393577204417721748144833));
            break;

          case 16:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(-0.09501250983763744018531933542495806313035305568));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(+0.09501250983763744018531933542495806313035305568));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(-0.28160355077925891323046050146049610648606949077));
            rule.get_coord ( 3) = Coord_(FEAT_F128C(+0.28160355077925891323046050146049610648606949077));
            rule.get_coord ( 4) = Coord_(FEAT_F128C(-0.45801677765722738634241944298357757354003161303));
            rule.get_coord ( 5) = Coord_(FEAT_F128C(+0.45801677765722738634241944298357757354003161303));
            rule.get_coord ( 6) = Coord_(FEAT_F128C(-0.61787624440264374844667176404879101899188221776));
            rule.get_coord ( 7) = Coord_(FEAT_F128C(+0.61787624440264374844667176404879101899188221776));
            rule.get_coord ( 8) = Coord_(FEAT_F128C(-0.75540440835500303389510119484744226835381365645));
            rule.get_coord ( 9) = Coord_(FEAT_F128C(+0.75540440835500303389510119484744226835381365645));
            rule.get_coord (10) = Coord_(FEAT_F128C(-0.86563120238783174388046789771239313238733538484));
            rule.get_coord (11) = Coord_(FEAT_F128C(+0.86563120238783174388046789771239313238733538484));
            rule.get_coord (12) = Coord_(FEAT_F128C(-0.94457502307323257607798841553460834509113927259));
            rule.get_coord (13) = Coord_(FEAT_F128C(+0.94457502307323257607798841553460834509113927259));
            rule.get_coord (14) = Coord_(FEAT_F128C(-0.98940093499164993259615417345033262742627407165));
            rule.get_coord (15) = Coord_(FEAT_F128C(+0.98940093499164993259615417345033262742627407165));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.189450610455068496285396723208283105146908988395));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.189450610455068496285396723208283105146908988395));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.182603415044923588866763667969219939383556223654));
            rule.get_weight( 3) = Weight_(FEAT_F128C(0.182603415044923588866763667969219939383556223654));
            rule.get_weight( 4) = Weight_(FEAT_F128C(0.169156519395002538189312079030359962211639473416));
            rule.get_weight( 5) = Weight_(FEAT_F128C(0.169156519395002538189312079030359962211639473416));
            rule.get_weight( 6) = Weight_(FEAT_F128C(0.149595988816576732081501730547478548970491068207));
            rule.get_weight( 7) = Weight_(FEAT_F128C(0.149595988816576732081501730547478548970491068207));
            rule.get_weight( 8) = Weight_(FEAT_F128C(0.124628971255533872052476282192016420144886859222));
            rule.get_weight( 9) = Weight_(FEAT_F128C(0.124628971255533872052476282192016420144886859222));
            rule.get_weight(10) = Weight_(FEAT_F128C(0.095158511682492784809925107602246226355263503183));
            rule.get_weight(11) = Weight_(FEAT_F128C(0.095158511682492784809925107602246226355263503183));
            rule.get_weight(12) = Weight_(FEAT_F128C(0.062253523938647892862843836994377694274986508352));
            rule.get_weight(13) = Weight_(FEAT_F128C(0.062253523938647892862843836994377694274986508352));
            rule.get_weight(14) = Weight_(FEAT_F128C(0.027152459411754094851780572456018103512267375566));
            rule.get_weight(15) = Weight_(FEAT_F128C(0.027152459411754094851780572456018103512267375566));
            break;

          case 17:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(+0.0));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(-0.17848418149584785585067749365406555747541933269));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(+0.17848418149584785585067749365406555747541933269));
            rule.get_coord ( 3) = Coord_(FEAT_F128C(-0.35123176345387631529718551709534600504053975157));
            rule.get_coord ( 4) = Coord_(FEAT_F128C(+0.35123176345387631529718551709534600504053975157));
            rule.get_coord ( 5) = Coord_(FEAT_F128C(-0.51269053708647696788624656862955187458292372241));
            rule.get_coord ( 6) = Coord_(FEAT_F128C(+0.51269053708647696788624656862955187458292372241));
            rule.get_coord ( 7) = Coord_(FEAT_F128C(-0.65767115921669076585030221664300233514780589147));
            rule.get_coord ( 8) = Coord_(FEAT_F128C(+0.65767115921669076585030221664300233514780589147));
            rule.get_coord ( 9) = Coord_(FEAT_F128C(-0.78151400389680140692523005552047605022397247274));
            rule.get_coord (10) = Coord_(FEAT_F128C(+0.78151400389680140692523005552047605022397247274));
            rule.get_coord (11) = Coord_(FEAT_F128C(-0.88023915372698590212295569448815569262341681793));
            rule.get_coord (12) = Coord_(FEAT_F128C(+0.88023915372698590212295569448815569262341681793));
            rule.get_coord (13) = Coord_(FEAT_F128C(-0.95067552176876776122271695789580302144338504655));
            rule.get_coord (14) = Coord_(FEAT_F128C(+0.95067552176876776122271695789580302144338504655));
            rule.get_coord (15) = Coord_(FEAT_F128C(-0.99057547531441733567543401994066527650778985045));
            rule.get_coord (16) = Coord_(FEAT_F128C(+0.99057547531441733567543401994066527650778985045));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.179446470356206525458265644261885621448780319897));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.176562705366992646325270990113197239150924418000));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.176562705366992646325270990113197239150924418000));
            rule.get_weight( 3) = Weight_(FEAT_F128C(0.168004102156450044509970663788323155021198128965));
            rule.get_weight( 4) = Weight_(FEAT_F128C(0.168004102156450044509970663788323155021198128965));
            rule.get_weight( 5) = Weight_(FEAT_F128C(0.154045761076810288081431594801958611940483058471));
            rule.get_weight( 6) = Weight_(FEAT_F128C(0.154045761076810288081431594801958611940483058471));
            rule.get_weight( 7) = Weight_(FEAT_F128C(0.135136368468525473286319981702350197372125853234));
            rule.get_weight( 8) = Weight_(FEAT_F128C(0.135136368468525473286319981702350197372125853234));
            rule.get_weight( 9) = Weight_(FEAT_F128C(0.111883847193403971094788385626355926735843424263));
            rule.get_weight(10) = Weight_(FEAT_F128C(0.111883847193403971094788385626355926735843424263));
            rule.get_weight(11) = Weight_(FEAT_F128C(0.085036148317179180883535370191062073850491389218));
            rule.get_weight(12) = Weight_(FEAT_F128C(0.085036148317179180883535370191062073850491389218));
            rule.get_weight(13) = Weight_(FEAT_F128C(0.055459529373987201129440165358244660512846251953));
            rule.get_weight(14) = Weight_(FEAT_F128C(0.055459529373987201129440165358244660512846251953));
            rule.get_weight(15) = Weight_(FEAT_F128C(0.024148302868547931960110026287565324691697315945));
            rule.get_weight(16) = Weight_(FEAT_F128C(0.024148302868547931960110026287565324691697315945));
            break;

          case 18:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(-0.08477501304173530124226185293578381173331738690));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(+0.08477501304173530124226185293578381173331738690));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(-0.25188622569150550958897285487791123016286176565));
            rule.get_coord ( 3) = Coord_(FEAT_F128C(+0.25188622569150550958897285487791123016286176565));
            rule.get_coord ( 4) = Coord_(FEAT_F128C(-0.41175116146284264603593179383305163707898968212));
            rule.get_coord ( 5) = Coord_(FEAT_F128C(+0.41175116146284264603593179383305163707898968212));
            rule.get_coord ( 6) = Coord_(FEAT_F128C(-0.55977083107394753460787154852532913692762648577));
            rule.get_coord ( 7) = Coord_(FEAT_F128C(+0.55977083107394753460787154852532913692762648577));
            rule.get_coord ( 8) = Coord_(FEAT_F128C(-0.69168704306035320787489108128884838945227057281));
            rule.get_coord ( 9) = Coord_(FEAT_F128C(+0.69168704306035320787489108128884838945227057281));
            rule.get_coord (10) = Coord_(FEAT_F128C(-0.80370495897252311568241745501459079710329892161));
            rule.get_coord (11) = Coord_(FEAT_F128C(+0.80370495897252311568241745501459079710329892161));
            rule.get_coord (12) = Coord_(FEAT_F128C(-0.89260246649755573920606059112714551540789527135));
            rule.get_coord (13) = Coord_(FEAT_F128C(+0.89260246649755573920606059112714551540789527135));
            rule.get_coord (14) = Coord_(FEAT_F128C(-0.95582394957139775518119589292977630997284413481));
            rule.get_coord (15) = Coord_(FEAT_F128C(+0.95582394957139775518119589292977630997284413481));
            rule.get_coord (16) = Coord_(FEAT_F128C(-0.99156516842093094673001600470615077025257893684));
            rule.get_coord (17) = Coord_(FEAT_F128C(+0.99156516842093094673001600470615077025257893684));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.169142382963143591840656470134986610334105819370));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.169142382963143591840656470134986610334105819370));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.164276483745832722986053776465927590412338953997));
            rule.get_weight( 3) = Weight_(FEAT_F128C(0.164276483745832722986053776465927590412338953997));
            rule.get_weight( 4) = Weight_(FEAT_F128C(0.154684675126265244925418003836374772193218396267));
            rule.get_weight( 5) = Weight_(FEAT_F128C(0.154684675126265244925418003836374772193218396267));
            rule.get_weight( 6) = Weight_(FEAT_F128C(0.140642914670650651204731303751947228095502410330));
            rule.get_weight( 7) = Weight_(FEAT_F128C(0.140642914670650651204731303751947228095502410330));
            rule.get_weight( 8) = Weight_(FEAT_F128C(0.122555206711478460184519126800201555228163897333));
            rule.get_weight( 9) = Weight_(FEAT_F128C(0.122555206711478460184519126800201555228163897333));
            rule.get_weight(10) = Weight_(FEAT_F128C(0.100942044106287165562813984924834607062801138887));
            rule.get_weight(11) = Weight_(FEAT_F128C(0.100942044106287165562813984924834607062801138887));
            rule.get_weight(12) = Weight_(FEAT_F128C(0.076425730254889056529129677616636525605317906208));
            rule.get_weight(13) = Weight_(FEAT_F128C(0.076425730254889056529129677616636525605317906208));
            rule.get_weight(14) = Weight_(FEAT_F128C(0.049714548894969796453334946202638641680866246128));
            rule.get_weight(15) = Weight_(FEAT_F128C(0.049714548894969796453334946202638641680866246128));
            rule.get_weight(16) = Weight_(FEAT_F128C(0.021616013526483310313342710266452469387685231475));
            rule.get_weight(17) = Weight_(FEAT_F128C(0.021616013526483310313342710266452469387685231475));
            break;

          case 19:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(+0.0));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(-0.16035864564022537586809611574074354950487350047));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(+0.16035864564022537586809611574074354950487350047));
            rule.get_coord ( 3) = Coord_(FEAT_F128C(-0.31656409996362983199011732884984491789228521913));
            rule.get_coord ( 4) = Coord_(FEAT_F128C(+0.31656409996362983199011732884984491789228521913));
            rule.get_coord ( 5) = Coord_(FEAT_F128C(-0.46457074137596094571726714810410236797628571462));
            rule.get_coord ( 6) = Coord_(FEAT_F128C(+0.46457074137596094571726714810410236797628571462));
            rule.get_coord ( 7) = Coord_(FEAT_F128C(-0.60054530466168102346963816494623927986832208273));
            rule.get_coord ( 8) = Coord_(FEAT_F128C(+0.60054530466168102346963816494623927986832208273));
            rule.get_coord ( 9) = Coord_(FEAT_F128C(-0.72096617733522937861709586082378162965714183290));
            rule.get_coord (10) = Coord_(FEAT_F128C(+0.72096617733522937861709586082378162965714183290));
            rule.get_coord (11) = Coord_(FEAT_F128C(-0.82271465653714282497892248671271390177453848620));
            rule.get_coord (12) = Coord_(FEAT_F128C(+0.82271465653714282497892248671271390177453848620));
            rule.get_coord (13) = Coord_(FEAT_F128C(-0.90315590361481790164266092853231248780939393405));
            rule.get_coord (14) = Coord_(FEAT_F128C(+0.90315590361481790164266092853231248780939393405));
            rule.get_coord (15) = Coord_(FEAT_F128C(-0.96020815213483003085277884068765152661509150327));
            rule.get_coord (16) = Coord_(FEAT_F128C(+0.96020815213483003085277884068765152661509150327));
            rule.get_coord (17) = Coord_(FEAT_F128C(-0.99240684384358440318901767025326049358931640140));
            rule.get_coord (18) = Coord_(FEAT_F128C(+0.99240684384358440318901767025326049358931640140));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.161054449848783695979163625320916735039902558578));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.158968843393954347649956439465047201678780158195));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.158968843393954347649956439465047201678780158195));
            rule.get_weight( 3) = Weight_(FEAT_F128C(0.152766042065859666778855400897662998461008267236));
            rule.get_weight( 4) = Weight_(FEAT_F128C(0.152766042065859666778855400897662998461008267236));
            rule.get_weight( 5) = Weight_(FEAT_F128C(0.142606702173606611775746109441902972475668344824));
            rule.get_weight( 6) = Weight_(FEAT_F128C(0.142606702173606611775746109441902972475668344824));
            rule.get_weight( 7) = Weight_(FEAT_F128C(0.128753962539336227675515784856877117055839577093));
            rule.get_weight( 8) = Weight_(FEAT_F128C(0.128753962539336227675515784856877117055839577093));
            rule.get_weight( 9) = Weight_(FEAT_F128C(0.111566645547333994716023901681765997481331853839));
            rule.get_weight(10) = Weight_(FEAT_F128C(0.111566645547333994716023901681765997481331853839));
            rule.get_weight(11) = Weight_(FEAT_F128C(0.091490021622449999464462094123839652660911651296));
            rule.get_weight(12) = Weight_(FEAT_F128C(0.091490021622449999464462094123839652660911651296));
            rule.get_weight(13) = Weight_(FEAT_F128C(0.069044542737641226580708258006013044961848031687));
            rule.get_weight(14) = Weight_(FEAT_F128C(0.069044542737641226580708258006013044961848031687));
            rule.get_weight(15) = Weight_(FEAT_F128C(0.044814226765699600332838157401994211951754227467));
            rule.get_weight(16) = Weight_(FEAT_F128C(0.044814226765699600332838157401994211951754227467));
            rule.get_weight(17) = Weight_(FEAT_F128C(0.019461788229726477036312041464438435752906609069));
            rule.get_weight(18) = Weight_(FEAT_F128C(0.019461788229726477036312041464438435752906609069));
            break;

          case 20:
            rule.get_coord ( 0) = Coord_(FEAT_F128C(-0.07652652113349733375464040939883821100479626681));
            rule.get_coord ( 1) = Coord_(FEAT_F128C(+0.07652652113349733375464040939883821100479626681));
            rule.get_coord ( 2) = Coord_(FEAT_F128C(-0.22778585114164507808049619536857462474308893768));
            rule.get_coord ( 3) = Coord_(FEAT_F128C(+0.22778585114164507808049619536857462474308893768));
            rule.get_coord ( 4) = Coord_(FEAT_F128C(-0.37370608871541956067254817702492723739574632170));
            rule.get_coord ( 5) = Coord_(FEAT_F128C(+0.37370608871541956067254817702492723739574632170));
            rule.get_coord ( 6) = Coord_(FEAT_F128C(-0.51086700195082709800436405095525099842549132920));
            rule.get_coord ( 7) = Coord_(FEAT_F128C(+0.51086700195082709800436405095525099842549132920));
            rule.get_coord ( 8) = Coord_(FEAT_F128C(-0.63605368072651502545283669622628593674338911679));
            rule.get_coord ( 9) = Coord_(FEAT_F128C(+0.63605368072651502545283669622628593674338911679));
            rule.get_coord (10) = Coord_(FEAT_F128C(-0.74633190646015079261430507035564159031073067956));
            rule.get_coord (11) = Coord_(FEAT_F128C(+0.74633190646015079261430507035564159031073067956));
            rule.get_coord (12) = Coord_(FEAT_F128C(-0.83911697182221882339452906170152068532962936506));
            rule.get_coord (13) = Coord_(FEAT_F128C(+0.83911697182221882339452906170152068532962936506));
            rule.get_coord (14) = Coord_(FEAT_F128C(-0.91223442825132590586775244120329811304918479742));
            rule.get_coord (15) = Coord_(FEAT_F128C(+0.91223442825132590586775244120329811304918479742));
            rule.get_coord (16) = Coord_(FEAT_F128C(-0.96397192727791379126766613119727722191206032780));
            rule.get_coord (17) = Coord_(FEAT_F128C(+0.96397192727791379126766613119727722191206032780));
            rule.get_coord (18) = Coord_(FEAT_F128C(-0.99312859918509492478612238847132027822264713090));
            rule.get_coord (19) = Coord_(FEAT_F128C(+0.99312859918509492478612238847132027822264713090));
            rule.get_weight( 0) = Weight_(FEAT_F128C(0.152753387130725850698084331955097593491948645112));
            rule.get_weight( 1) = Weight_(FEAT_F128C(0.152753387130725850698084331955097593491948645112));
            rule.get_weight( 2) = Weight_(FEAT_F128C(0.149172986472603746787828737001969436692679904081));
            rule.get_weight( 3) = Weight_(FEAT_F128C(0.149172986472603746787828737001969436692679904081));
            rule.get_weight( 4) = Weight_(FEAT_F128C(0.142096109318382051329298325067164933034515413392));
            rule.get_weight( 5) = Weight_(FEAT_F128C(0.142096109318382051329298325067164933034515413392));
            rule.get_weight( 6) = Weight_(FEAT_F128C(0.131688638449176626898494499748163134916110511146));
            rule.get_weight( 7) = Weight_(FEAT_F128C(0.131688638449176626898494499748163134916110511146));
            rule.get_weight( 8) = Weight_(FEAT_F128C(0.118194531961518417312377377711382287005041219548));
            rule.get_weight( 9) = Weight_(FEAT_F128C(0.118194531961518417312377377711382287005041219548));
            rule.get_weight(10) = Weight_(FEAT_F128C(0.101930119817240435036750135480349876166691656023));
            rule.get_weight(11) = Weight_(FEAT_F128C(0.101930119817240435036750135480349876166691656023));
            rule.get_weight(12) = Weight_(FEAT_F128C(0.083276741576704748724758143222046206100177828583));
            rule.get_weight(13) = Weight_(FEAT_F128C(0.083276741576704748724758143222046206100177828583));
            rule.get_weight(14) = Weight_(FEAT_F128C(0.062672048334109063569506535187041606351601076578));
            rule.get_weight(15) = Weight_(FEAT_F128C(0.062672048334109063569506535187041606351601076578));
            rule.get_weight(16) = Weight_(FEAT_F128C(0.040601429800386941331039952274932109879090639989));
            rule.get_weight(17) = Weight_(FEAT_F128C(0.040601429800386941331039952274932109879090639989));
            rule.get_weight(18) = Weight_(FEAT_F128C(0.017614007139152118311861962351852816362143105543));
            rule.get_weight(19) = Weight_(FEAT_F128C(0.017614007139152118311861962351852816362143105543));
            break;
          }
        }
      }; // class GaussLegendreDriver<...>
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_SCALAR_GAUSS_LEGENDRE_DRIVER_HPP
