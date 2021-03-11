// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_CUBATURE_DUNAVANT_DRIVER_HPP
#define KERNEL_CUBATURE_DUNAVANT_DRIVER_HPP 1

// includes, FEAT
#include <kernel/cubature/symmetric_simplex_driver.hpp>
#include <kernel/util/meta_math.hpp>

namespace FEAT
{
  namespace Cubature
  {
    /**
     * \brief Dunavant "open" driver class template
     *
     * This driver implements the open Dunavant rules for triangles.
     * \see D.A. Dunavant: High Degree Efficient Symmetrical Gaussian Quadrature Rules for the Triangle
     *      Int. j. numer. methods eng., Volume 21 (1985), pp. 1129 - 1148
     *
     * \tparam Shape_
     * The shape type of the element.
     *
     * \author Peter Zajac, Stefan Wahlers
     */
    template<typename Shape_>
    class DunavantDriver DOXY({});

    // Simplex<2> specialization
    template<>
    class DunavantDriver<Shape::Simplex<2> > :
      public SymmetricSimplexDriver<Shape::Simplex<2> >
    {
    public:
      /// this rule is variadic
      static constexpr bool variadic = true;
      static constexpr int min_points = 2;
      static constexpr int max_points = 20;

      /// Returns the name of the cubature rule.
      static String name()
      {
        return "dunavant";
      }

      static int count(int points)
      {
        switch(points)
        {
        case 2:
          return 3;
        case 3:
          return 4;
        case 4:
          return 6;
        case 5:
          return 7;
        case 6:
          return 12;
        case 7:
          return 13;
        case 8:
          return 16;
        case 9:
          return 19;
        case 10:
          return 25;
        case 11:
          return 27;
        case 12:
          return 33;
        case 13:
          return 37;
        case 14:
          return 42;
        case 15:
          return 48;
        case 16:
          return 52;
        case 17:
          return 61;
        case 18:
          return 70;
        case 19:
          return 73;
        case 20:
          return 79;
        default:
          return 0;
        }
      }

      /**
       * \brief Fills the cubature rule structure.
       *
       * \param[in,out] rule
       * The cubature rule to be filled.
       */
      template<
        typename Weight_,
        typename Coord_,
        typename Point_>
      static void fill(Rule<Shape::Simplex<2>, Weight_, Coord_, Point_>& rule, int num_points)
      {
        int off(0);
        switch(num_points)
        {
        case 2:
          off += fill_sym2(rule, off, Weight_(0.1666666666666665), Coord_(0.666666666666667), Coord_(0.166666666666667));
          break;
        case 3:
          off += fill_sym1(rule, off, Weight_(-0.2812500000000000), Coord_(0.333333333333333));
          off += fill_sym2(rule, off, Weight_(0.2604166666666665), Coord_(0.600000000000000), Coord_(0.200000000000000));
          break;
        case 4:
          off += fill_sym2(rule, off, Weight_(0.1116907948390055), Coord_(0.108103018168070), Coord_(0.445948490915965));
          off += fill_sym2(rule, off, Weight_(0.0549758718276610), Coord_(0.816847572980459), Coord_(0.091576213509771));
          break;
        case 5:
          off += fill_sym1(rule, off, Weight_(0.1125000000000000), Coord_(0.333333333333333));
          off += fill_sym2(rule, off, Weight_(0.0661970763942530), Coord_(0.059715871789770), Coord_(0.470142064105115));
          off += fill_sym2(rule, off, Weight_(0.0629695902724135), Coord_(0.797426985353087), Coord_(0.101286507323456));
          break;
        case 6:
          off += fill_sym2(rule, off, Weight_(0.0583931378631895), Coord_(0.501426509658179), Coord_(0.249286745170910));
          off += fill_sym2(rule, off, Weight_(0.0254224531851035), Coord_(0.873821971016996), Coord_(0.063089014491502));
          off += fill_sym3(rule, off, Weight_(0.0414255378091870), Coord_(0.053145049844817), Coord_(0.310352451033784), Coord_(0.636502499121399));
          break;
        case 7:
          off += fill_sym1(rule, off, Weight_(-0.149570044467682), Coord_(0.333333333333333));
          off += fill_sym2(rule, off, Weight_(0.0878076287166040), Coord_(0.479308067841920), Coord_(0.260345966079040));
          off += fill_sym2(rule, off, Weight_(0.0266736178044190), Coord_(0.869739794195568), Coord_(0.065130102902216));
          off += fill_sym3(rule, off, Weight_(0.0385568804451285), Coord_(0.048690315425316), Coord_(0.312865496004874), Coord_(0.638444188569810));
          break;
        case 8:
          off += fill_sym1(rule, off, Weight_(0.0721578038388935), Coord_(0.333333333333333));
          off += fill_sym2(rule, off, Weight_(0.0475458171336425), Coord_(0.081414823414554), Coord_(0.459292588292723));
          off += fill_sym2(rule, off, Weight_(0.0516086852673590), Coord_(0.658861384496480), Coord_(0.170569307751760));
          off += fill_sym2(rule, off, Weight_(0.0162292488115990), Coord_(0.898905543365938), Coord_(0.050547228317031));
          off += fill_sym3(rule, off, Weight_(0.0136151570872175), Coord_(0.008394777409958), Coord_(0.263112829634638), Coord_(0.728492392955404));
          break;
        case 9:
          off += fill_sym1(rule, off, Weight_(0.0485678981413995), Coord_(0.333333333333333));
          off += fill_sym2(rule, off, Weight_(0.0156673501135695), Coord_(0.020634961602525), Coord_(0.489682519198738));
          off += fill_sym2(rule, off, Weight_(0.0389137705023870), Coord_(0.125820817014127), Coord_(0.437089591492937));
          off += fill_sym2(rule, off, Weight_(0.0398238694636050), Coord_(0.623592928761935), Coord_(0.188203535619033));
          off += fill_sym2(rule, off, Weight_(0.0127888378293490), Coord_(0.910540973211095), Coord_(0.044729513394453));
          off += fill_sym3(rule, off, Weight_(0.0216417696886445), Coord_(0.036838412054736), Coord_(0.221962989160766), Coord_(0.741198598784498));
          break;
        case 10:
          off += fill_sym1(rule, off, Weight_(0.0454089951913770), Coord_(0.333333333333333));
          off += fill_sym2(rule, off, Weight_(0.0183629788782335), Coord_(0.028844733232685), Coord_(0.485577633383657));
          off += fill_sym2(rule, off, Weight_(0.0226605297177640), Coord_(0.781036849029926), Coord_(0.109481575485037));
          off += fill_sym3(rule, off, Weight_(0.0363789584227100), Coord_(0.141707219414880), Coord_(0.307939838764121), Coord_(0.550352941820999));
          off += fill_sym3(rule, off, Weight_(0.0141636212655285), Coord_(0.025003534762686), Coord_(0.246672560639903), Coord_(0.728323904597411));
          off += fill_sym3(rule, off, Weight_(0.0047108334818665), Coord_(0.009540815400299), Coord_(0.066803251012200), Coord_(0.923655933587500));
          break;
        case 11:
          off += fill_sym2(rule, off, Weight_(0.0004635031644805), Coord_(-0.069222096541517), Coord_(0.534611048270758));
          off += fill_sym2(rule, off, Weight_(0.0385747674574065), Coord_(0.202061394068290), Coord_(0.398969302965855));
          off += fill_sym2(rule, off, Weight_(0.0296614886903870), Coord_(0.593380199137435), Coord_(0.203309900431282));
          off += fill_sym2(rule, off, Weight_(0.0180922702517090), Coord_(0.761298175434837), Coord_(0.119350912282581));
          off += fill_sym2(rule, off, Weight_(0.0068298655013390), Coord_(0.935270103777448), Coord_(0.032364948111276));
          off += fill_sym3(rule, off, Weight_(0.0261685559811020), Coord_(0.050178138310495), Coord_(0.356620648261293), Coord_(0.593201213428213));
          off += fill_sym3(rule, off, Weight_(0.0103538298195705), Coord_(0.021022016536166), Coord_(0.171488980304042), Coord_(0.807489003159792));
          break;
        case 12:
          off += fill_sym2(rule, off, Weight_(0.0128655332202275), Coord_(0.023565220452390), Coord_(0.488217389773805));
          off += fill_sym2(rule, off, Weight_(0.0218462722690190), Coord_(0.120551215411079), Coord_(0.439724392294460));
          off += fill_sym2(rule, off, Weight_(0.0314291121089425), Coord_(0.457579229975768), Coord_(0.271210385012116));
          off += fill_sym2(rule, off, Weight_(0.0173980564653545), Coord_(0.744847708916828), Coord_(0.127576145541586));
          off += fill_sym2(rule, off, Weight_(0.0030831305257795), Coord_(0.957365299093579), Coord_(0.021317350453210));
          off += fill_sym3(rule, off, Weight_(0.0201857788831905), Coord_(0.115343494534698), Coord_(0.275713269685514), Coord_(0.608943235779788));
          off += fill_sym3(rule, off, Weight_(0.0111783866011515), Coord_(0.022838332222257), Coord_(0.281325580989940), Coord_(0.695836086787803));
          off += fill_sym3(rule, off, Weight_(0.0086581155543295), Coord_(0.025734050548330), Coord_(0.116251915907597), Coord_(0.858014033544073));
          break;
        case 13:
          off += fill_sym1(rule, off, Weight_(0.0262604617004010), Coord_(0.333333333333333));
          off += fill_sym2(rule, off, Weight_(0.0056400726046650), Coord_(0.009903630120591), Coord_(0.495048184939705));
          off += fill_sym2(rule, off, Weight_(0.0157117591812270), Coord_(0.062566729780852), Coord_(0.468716635109574));
          off += fill_sym2(rule, off, Weight_(0.0235362512520970), Coord_(0.170957326397447), Coord_(0.414521336801277));
          off += fill_sym2(rule, off, Weight_(0.0236817932681775), Coord_(0.541200855914337), Coord_(0.229399572042831));
          off += fill_sym2(rule, off, Weight_(0.0155837645228970), Coord_(0.771151009607340), Coord_(0.114424495196330));
          off += fill_sym2(rule, off, Weight_(0.0039878857325370), Coord_(0.950377217273082), Coord_(0.024811391363459));
          off += fill_sym3(rule, off, Weight_(0.0184242013643660), Coord_(0.094853828379579), Coord_(0.268794997058761), Coord_(0.636351174561660));
          off += fill_sym3(rule, off, Weight_(0.0087007316519110), Coord_(0.018100773278807), Coord_(0.291730066734288), Coord_(0.690169159986905));
          off += fill_sym3(rule, off, Weight_(0.0077608934195225), Coord_(0.022233076674090), Coord_(0.126357385491669), Coord_(0.85140953783424));
          break;
        case 14:
          off += fill_sym2(rule, off, Weight_(0.0109417906847145), Coord_(0.022072179275643), Coord_(0.488963910362179));
          off += fill_sym2(rule, off, Weight_(0.0163941767720625), Coord_(0.164710561319092), Coord_(0.417644719340454));
          off += fill_sym2(rule, off, Weight_(0.0258870522536460), Coord_(0.453044943382323), Coord_(0.273477528308839));
          off += fill_sym2(rule, off, Weight_(0.0210812943684965), Coord_(0.645588935174913), Coord_(0.177205532412543));
          off += fill_sym2(rule, off, Weight_(0.0072168498348885), Coord_(0.876400233818255), Coord_(0.061799883090873));
          off += fill_sym2(rule, off, Weight_(0.0024617018012000), Coord_(0.961218077502598), Coord_(0.019390961248701));
          off += fill_sym3(rule, off, Weight_(0.0123328766062820), Coord_(0.057124757403648), Coord_(0.172266687821356), Coord_(0.770608554774996));
          off += fill_sym3(rule, off, Weight_(0.0192857553935305), Coord_(0.092916249356972), Coord_(0.336861459796345), Coord_(0.570222290846683));
          off += fill_sym3(rule, off, Weight_(0.0072181540567670), Coord_(0.014646950055654), Coord_(0.298372882136258), Coord_(0.686980167808088));
          off += fill_sym3(rule, off, Weight_(0.0025051144192505), Coord_(0.001268330932872), Coord_(0.118974497696957), Coord_(0.879757171370171));
          break;
        case 15:
          off += fill_sym2(rule, off, Weight_(0.0009584378214245), Coord_(-0.013945833716486), Coord_(0.506972916858243));
          off += fill_sym2(rule, off, Weight_(0.0221245136355725), Coord_(0.137187291433955), Coord_(0.431406354283023));
          off += fill_sym2(rule, off, Weight_(0.0255932743594260), Coord_(0.444612710305711), Coord_(0.277693644847144));
          off += fill_sym2(rule, off, Weight_(0.0118438679353440), Coord_(0.747070217917492), Coord_(0.126464891041254));
          off += fill_sym2(rule, off, Weight_(0.0066448878450105), Coord_(0.858383228050628), Coord_(0.070808385974686));
          off += fill_sym2(rule, off, Weight_(0.0023744583040960), Coord_(0.962069659517853), Coord_(0.018965170241073));
          off += fill_sym3(rule, off, Weight_(0.0192750362997965), Coord_(0.133734161966621), Coord_(0.261311371140087), Coord_(0.604954466893291));
          off += fill_sym3(rule, off, Weight_(0.0136079071603120), Coord_(0.036366677396917), Coord_(0.388046767090269), Coord_(0.575586555512814));
          off += fill_sym3(rule, off, Weight_(0.0010910386833985), Coord_(-0.010174883126571), Coord_(0.285712220049916), Coord_(0.724462663076655));
          off += fill_sym3(rule, off, Weight_(0.0107526599238655), Coord_(0.036843869875878), Coord_(0.215599664072284), Coord_(0.747556466051838));
          off += fill_sym3(rule, off, Weight_(0.0038369713155245), Coord_(0.012459809331199), Coord_(0.103575616576386), Coord_(0.883964574092416));
          break;
        case 16:
          off += fill_sym1(rule, off, Weight_(0.0234378487138210), Coord_(0.333333333333333));
          off += fill_sym2(rule, off, Weight_(0.0032029392892925), Coord_(0.005238916103123), Coord_(0.497380541948438));
          off += fill_sym2(rule, off, Weight_(0.0208551483696935), Coord_(0.173061122901295), Coord_(0.413469438549352));
          off += fill_sym2(rule, off, Weight_(0.0134457421250320), Coord_(0.059082801866017), Coord_(0.470458599066991));
          off += fill_sym2(rule, off, Weight_(0.0210662613808250), Coord_(0.518892500060958), Coord_(0.240553749969521));
          off += fill_sym2(rule, off, Weight_(0.0150001334213865), Coord_(0.704068411554854), Coord_(0.147965794222573));
          off += fill_sym2(rule, off, Weight_(0.0071000494625120), Coord_(0.849069624685052), Coord_(0.075465187657474));
          off += fill_sym2(rule, off, Weight_(0.0017912311756365), Coord_(0.966807194753950), Coord_(0.016596402623025));
          off += fill_sym3(rule, off, Weight_(0.0163865737303135), Coord_(0.103575692245252), Coord_(0.296555596579887), Coord_(0.599868711174861));
          off += fill_sym3(rule, off, Weight_(0.0076491531242205), Coord_(0.020083411655416), Coord_(0.337723063403079), Coord_(0.642193524941505));
          off += fill_sym3(rule, off, Weight_(0.0011931220964195), Coord_(-0.004341002614139), Coord_(0.204748281642812), Coord_(0.799592720971327));
          off += fill_sym3(rule, off, Weight_(0.0095423963779495), Coord_(0.041941786468010), Coord_(0.189358492130623), Coord_(0.768699721401368));
          off += fill_sym3(rule, off, Weight_(0.0034250272732710), Coord_(0.014317320230681), Coord_(0.085283615682657), Coord_(0.900399064086661));
          break;
        case 17:
          off += fill_sym1(rule, off, Weight_(0.0167185996454015), Coord_(0.333333333333333));
          off += fill_sym2(rule, off, Weight_(0.0025467077202535), Coord_(0.005658918886452), Coord_(0.497170540556774));
          off += fill_sym2(rule, off, Weight_(0.0073354322638190), Coord_(0.035647354750751), Coord_(0.482176322624625));
          off += fill_sym2(rule, off, Weight_(0.0121754391768360), Coord_(0.099520061958437), Coord_(0.450239969020782));
          off += fill_sym2(rule, off, Weight_(0.0155537754344845), Coord_(0.199467521245206), Coord_(0.400266239377397));
          off += fill_sym2(rule, off, Weight_(0.0156285556093100), Coord_(0.495717464058095), Coord_(0.252141267970953));
          off += fill_sym2(rule, off, Weight_(0.0124078271698325), Coord_(0.675905990683077), Coord_(0.162047004658461));
          off += fill_sym2(rule, off, Weight_(0.0070280365352785), Coord_(0.848248235478508), Coord_(0.075875882260746));
          off += fill_sym2(rule, off, Weight_(0.0015973380868895), Coord_(0.968690546064356), Coord_(0.015654726967822));
          off += fill_sym3(rule, off, Weight_(0.0040598276594965), Coord_(0.010186928826919), Coord_(0.334319867363658), Coord_(0.655493203809423));
          off += fill_sym3(rule, off, Weight_(0.0134028711415815), Coord_(0.135440871671036), Coord_(0.292221537796944), Coord_(0.572337590532020));
          off += fill_sym3(rule, off, Weight_(0.0092299966054110), Coord_(0.054423924290583), Coord_(0.319574885423190), Coord_(0.626001190286228));
          off += fill_sym3(rule, off, Weight_(0.0042384342671640), Coord_(0.012868560833637), Coord_(0.190704224192292), Coord_(0.796427214974071));
          off += fill_sym3(rule, off, Weight_(0.0091463983850125), Coord_(0.067165782413524), Coord_(0.180483211648746), Coord_(0.752351005937729));
          off += fill_sym3(rule, off, Weight_(0.0033328160020825), Coord_(0.014663182224828), Coord_(0.080711313679564), Coord_(0.904625504095608));
          break;
        case 18:
          off += fill_sym1(rule, off, Weight_(0.0154049699688235), Coord_(0.333333333333333));
          off += fill_sym2(rule, off, Weight_(0.0045362183397020), Coord_(0.013310382738157), Coord_(0.493344808630921));
          off += fill_sym2(rule, off, Weight_(0.0093806584697970), Coord_(0.061578811516086), Coord_(0.469210594241957));
          off += fill_sym2(rule, off, Weight_(0.0097205489927385), Coord_(0.127437208225989), Coord_(0.436261395887006));
          off += fill_sym2(rule, off, Weight_(0.0138769743054050), Coord_(0.210307658653168), Coord_(0.394846170673416));
          off += fill_sym2(rule, off, Weight_(0.0161281126757285), Coord_(0.500410862393686), Coord_(0.249794568803157));
          off += fill_sym2(rule, off, Weight_(0.0125370163084610), Coord_(0.677135612512315), Coord_(0.161432193743843));
          off += fill_sym2(rule, off, Weight_(0.0076359639859160), Coord_(0.846803545029257), Coord_(0.076598227485371));
          off += fill_sym2(rule, off, Weight_(0.0033969610114815), Coord_(0.951495121293100), Coord_(0.024252439353450));
          off += fill_sym2(rule, off, Weight_(-0.0011115493649600), Coord_(0.913707265566071), Coord_(0.043146367216965));
          off += fill_sym3(rule, off, Weight_(0.0031659570382030), Coord_(0.008430536202420), Coord_(0.358911494940944), Coord_(0.632657968856636));
          off += fill_sym3(rule, off, Weight_(0.0136287690245690), Coord_(0.131186551737188), Coord_(0.294402476751957), Coord_(0.574410971510855));
          off += fill_sym3(rule, off, Weight_(0.0088383928247325), Coord_(0.050203151565675), Coord_(0.325017801641814), Coord_(0.624779046792512));
          off += fill_sym3(rule, off, Weight_(0.0091897423190350), Coord_(0.066329263810916), Coord_(0.184737559666046), Coord_(0.748933176523037));
          off += fill_sym3(rule, off, Weight_(0.0040523664040960), Coord_(0.011996194566236), Coord_(0.218796800013321), Coord_(0.769207005420443));
          off += fill_sym3(rule, off, Weight_(0.0038170645353625), Coord_(0.014858100590125), Coord_(0.101179597136408), Coord_(0.883962302273467));
          off += fill_sym3(rule, off, Weight_(0.0000230938303970), Coord_(-0.035222015287949), Coord_(0.020874755282586), Coord_(1.014347260005363));
          break;
        case 19:
          off += fill_sym1(rule, off, Weight_(0.0164531656944595), Coord_(0.333333333333333));
          off += fill_sym2(rule, off, Weight_(0.0051653659456360), Coord_(0.020780025853987), Coord_(0.489609987073006));
          off += fill_sym2(rule, off, Weight_(0.0111936236315080), Coord_(0.090926214604215), Coord_(0.454536892697893));
          off += fill_sym2(rule, off, Weight_(0.0151330629347340), Coord_(0.197166638701138), Coord_(0.401416680649431));
          off += fill_sym2(rule, off, Weight_(0.0152454839010990), Coord_(0.488896691193805), Coord_(0.255551654403098));
          off += fill_sym2(rule, off, Weight_(0.0120796063708205), Coord_(0.645844115695741), Coord_(0.177077942152130));
          off += fill_sym2(rule, off, Weight_(0.0080254017934005), Coord_(0.779877893544096), Coord_(0.110061053227952));
          off += fill_sym2(rule, off, Weight_(0.0040422901308920), Coord_(0.888942751496321), Coord_(0.055528624251840));
          off += fill_sym2(rule, off, Weight_(0.0010396810137425), Coord_(0.974756272445543), Coord_(0.012621863777229));
          off += fill_sym3(rule, off, Weight_(0.0019424384524905), Coord_(0.003611417848412), Coord_(0.395754787356943), Coord_(0.600633794794645));
          off += fill_sym3(rule, off, Weight_(0.0127870803060110), Coord_(0.134466754530780), Coord_(0.307929983880436), Coord_(0.557603261588784));
          off += fill_sym3(rule, off, Weight_(0.0044404517866690), Coord_(0.014446025776115), Coord_(0.264566948406520), Coord_(0.720987025817365));
          off += fill_sym3(rule, off, Weight_(0.0080622733808655), Coord_(0.046933578838178), Coord_(0.358539352205951), Coord_(0.594527068955871));
          off += fill_sym3(rule, off, Weight_(0.0012459709087455), Coord_(0.002861120350567), Coord_(0.157807405968595), Coord_(0.839331473680839));
          off += fill_sym3(rule, off, Weight_(0.0091214200594755), Coord_(0.223861424097916), Coord_(0.075050596975911), Coord_(0.701087978926173));
          off += fill_sym3(rule, off, Weight_(0.0051292818680995), Coord_(0.034647074816760), Coord_(0.142421601113383), Coord_(0.822931324069857));
          off += fill_sym3(rule, off, Weight_(0.0018999644276510), Coord_(0.010161119296278), Coord_(0.065494628082938), Coord_(0.924344252620784));
          break;
        case 20:
          off += fill_sym1(rule, off, Weight_(0.0165285277708120), Coord_(0.333333333333333));
          off += fill_sym2(rule, off, Weight_(0.0004335095928315), Coord_(-0.001900928704400), Coord_(0.500950464352200));
          off += fill_sym2(rule, off, Weight_(0.0058300263582240), Coord_(0.023574084130543), Coord_(0.488212957934729));
          off += fill_sym2(rule, off, Weight_(0.0114384681782105), Coord_(0.089726636099435), Coord_(0.455136681950283));
          off += fill_sym2(rule, off, Weight_(0.0152244913369690), Coord_(0.196007481363421), Coord_(0.401996259318289));
          off += fill_sym2(rule, off, Weight_(0.0153124458626775), Coord_(0.488214180481157), Coord_(0.255892909759421));
          off += fill_sym2(rule, off, Weight_(0.0121840288384000), Coord_(0.647023488009788), Coord_(0.176488255995106));
          off += fill_sym2(rule, off, Weight_(0.0079987160160120), Coord_(0.791658289326483), Coord_(0.104170855336758));
          off += fill_sym2(rule, off, Weight_(0.0038491509078010), Coord_(0.893862072318140), Coord_(0.053068963840930));
          off += fill_sym2(rule, off, Weight_(-0.0003160302487440), Coord_(0.916762569607942), Coord_(0.041618715196029));
          off += fill_sym2(rule, off, Weight_(0.0008755671505965), Coord_(0.976836157186356), Coord_(0.011581921406822));
          off += fill_sym3(rule, off, Weight_(0.0082329195947880), Coord_(0.048741583664839), Coord_(0.344855770229001), Coord_(0.606402646106160));
          off += fill_sym3(rule, off, Weight_(0.0024195167702425), Coord_(0.006314115948605), Coord_(0.377843269594854), Coord_(0.615842614456541));
          off += fill_sym3(rule, off, Weight_(0.0129024532673250), Coord_(0.134316520547348), Coord_(0.306635479062357), Coord_(0.559048000390295));
          off += fill_sym3(rule, off, Weight_(0.0042355455272205), Coord_(0.013973893962392), Coord_(0.249419362774742), Coord_(0.736606743262866));
          off += fill_sym3(rule, off, Weight_(0.0091774570531400), Coord_(0.075549132909764), Coord_(0.212775724802802), Coord_(0.711675142287434));
          off += fill_sym3(rule, off, Weight_(0.0003522023389540), Coord_(-0.008368153208227), Coord_(0.146965436053239), Coord_(0.861402717154987));
          off += fill_sym3(rule, off, Weight_(0.0050563424637310), Coord_(0.026686063258714), Coord_(0.137726978828923), Coord_(0.835586957912363));
          off += fill_sym3(rule, off, Weight_(0.0017869546929750), Coord_(0.010547719294141), Coord_(0.059696109149007), Coord_(0.929756171556853));
          break;
        }
      }
    }; // class DunavantDriver<Simplex<...>,...>

  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_DUNAVANT_DRIVER_HPP
