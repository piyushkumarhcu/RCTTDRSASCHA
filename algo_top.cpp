#include "algo_top.h"
#include "algo_top_parameters.h"
#include "bitonicSort16.h"

void processInputLinks(ap_uint<576> link_in[N_INPUT_LINKS], crystal ECALRegion3x4_1[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], crystal ECALRegion3x4_2[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI],crystal ECALRegion3x4_3[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI]){
	#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
        #pragma HLS ARRAY_PARTITION variable=ECALRegion3x4_1 complete dim=0
        #pragma HLS ARRAY_PARTITION variable=ECALRegion3x4_2 complete dim=0
        #pragma HLS ARRAY_PARTITION variable=ECALRegion3x4_3 complete dim=0

        ap_uint<32> start = 0;
        ap_uint<32> end = 13;

        ap_uint<6> wordId, wordId1, wordId2, startId ;


        for(loop i=0; i<CRYSTAL_IN_ETA; i++){
//                 #pragma HLS UNROLL
        	for(loop j=0; j<CRYSTAL_IN_PHI; j++){
            	#pragma HLS UNROLL
        		wordId 	= (i/5)*4+(j/5);
                wordId1 = wordId + 12;
                wordId2 = wordId + 24;
                startId = (i%5)*5+(j%5) ;
                start 	= startId*14 ; end = start + 13 ;

                ECALRegion3x4_1[i][j] = crystal(link_in[wordId].range(end, start));
                ECALRegion3x4_2[i][j] = crystal(link_in[wordId1].range(end, start));
                if(wordId2 <= 31) ECALRegion3x4_3[i][j] = crystal(link_in[wordId2].range(end, start));
                else ECALRegion3x4_3[i][j] = crystal(0);
             }
        }
}

void processOutLink(Cluster inCluster[3], tower_t inTower[16], ap_uint<576> &outLink){
        #pragma HLS PIPELINE II=6
        #pragma HLS LATENCY min=6

        ap_uint<32> start = 0;

        outLink = 0;

        for(loop oLink=0; oLink<16; oLink++){
        	#pragma HLS UNROLL
                ap_uint<32> end = start + 11;
                outLink.range(end, start) = inTower[oLink].et();
                start += 12;
        }

        for(loop oLink=0; oLink<3; oLink++){
        	#pragma HLS UNROLL
        	ap_uint<32> end = start + 59;
        	outLink.range(end, start) = inCluster[oLink];
        	start += 60;
        }

}

ecaltp_t bestOf2(const ecaltp_t& ecaltp0, const ecaltp_t& ecaltp1) {
        ecaltp_t x;
        x = (ecaltp0.energy > ecaltp1.energy)?ecaltp0:ecaltp1;

        return x;
}

ecaltp_t getPeakBin20N(const etaStrip_t& etaStrip){
//#pragma HLS PIPELINE II=2
#pragma HLS latency min=4

ecaltp_t best01       = bestOf2(etaStrip.cr0,etaStrip.cr1) ;
ecaltp_t best23       = bestOf2(etaStrip.cr2,etaStrip.cr3) ;
ecaltp_t best45       = bestOf2(etaStrip.cr4,etaStrip.cr5) ;
ecaltp_t best67       = bestOf2(etaStrip.cr6,etaStrip.cr7) ;
ecaltp_t best89       = bestOf2(etaStrip.cr8,etaStrip.cr9) ;
ecaltp_t best1011     = bestOf2(etaStrip.cr10,etaStrip.cr11) ;
ecaltp_t best1213     = bestOf2(etaStrip.cr12,etaStrip.cr13) ;
ecaltp_t best1415     = bestOf2(etaStrip.cr14,etaStrip.cr15) ;
ecaltp_t best1617     = bestOf2(etaStrip.cr16,etaStrip.cr17) ;
ecaltp_t best1819     = bestOf2(etaStrip.cr18,etaStrip.cr19) ;

ecaltp_t best0123     = bestOf2(best01,best23) ;
ecaltp_t best4567     = bestOf2(best45,best67) ;
ecaltp_t best891011   = bestOf2(best89,best1011) ;
ecaltp_t best12131415 = bestOf2(best1213,best1415) ;
ecaltp_t best16171819 = bestOf2(best1617,best1819) ;

ecaltp_t best01234567 = bestOf2(best0123,best4567) ;
ecaltp_t best89101112131415 = bestOf2(best891011,best12131415) ;

ecaltp_t best0to15 		= bestOf2(best01234567,best89101112131415) ;
ecaltp_t bestOf20 		= bestOf2(best0to15,best16171819) ;

return bestOf20 ;
}

crystalMax getPeakBin15N(const etaStripPeak_t& etaStrip){
//#pragma HLS PIPELINE II=2
#pragma HLS latency min=4

crystalMax x;

ecaltp_t best01 = bestOf2(etaStrip.pk0,etaStrip.pk1) ;
ecaltp_t best23 = bestOf2(etaStrip.pk2,etaStrip.pk3) ;
ecaltp_t best45 = bestOf2(etaStrip.pk4,etaStrip.pk5) ;
ecaltp_t best67 = bestOf2(etaStrip.pk6,etaStrip.pk7) ;
ecaltp_t best89 = bestOf2(etaStrip.pk8,etaStrip.pk9) ;
ecaltp_t best1011 = bestOf2(etaStrip.pk10,etaStrip.pk11) ;
ecaltp_t best1213 = bestOf2(etaStrip.pk12,etaStrip.pk13) ;

ecaltp_t best0123 = bestOf2(best01,best23) ;
ecaltp_t best4567 = bestOf2(best45,best67) ;
ecaltp_t best891011 = bestOf2(best89,best1011) ;
ecaltp_t best121314 = bestOf2(best1213,etaStrip.pk14) ;

ecaltp_t best01234567 = bestOf2(best0123,best4567);
ecaltp_t best891011121314 = bestOf2(best891011,best121314) ;

ecaltp_t bestOf15 = bestOf2(best01234567,best891011121314) ;

        x.energy = bestOf15.energy ;
        x.etaMax = bestOf15.eta ;
        x.phiMax = bestOf15.phi ;

return x ;
}

ap_uint<12> getTowerEt(ap_uint<10> temp[5][5]){
//        #pragma HLS PIPELINE II=9
        #pragma HLS latency min=4
        #pragma HLS ARRAY_PARTITION variable=temp complete dim=0

		ap_uint<15> eta_strip[5];
		ap_uint<15> tEtSum;

		for(loop i=0; i<5 ; i++){
			#pragma HLS unroll
			eta_strip[i] = temp[i][0] + temp[i][1] + temp[i][2] + temp[i][3] + temp[i][4];
		}

		tEtSum = eta_strip[0] + eta_strip[1] + eta_strip[2] + eta_strip[3] + eta_strip[4];

		ap_uint<12> towerEt   = (tEtSum > 0xFFF) ? (ap_uint<12>)0xFFF : (ap_uint<12>) tEtSum;

		return towerEt;
}


clusterInfo getClusterPosition(const ecalRegion_t& ecalRegion){
 #pragma HLS latency min=4

        etaStripPeak_t etaStripPeak;
        clusterInfo cluster ;


                etaStripPeak.pk0  = getPeakBin20N(ecalRegion.etaStrip0);
                etaStripPeak.pk1  = getPeakBin20N(ecalRegion.etaStrip1);
                etaStripPeak.pk2  = getPeakBin20N(ecalRegion.etaStrip2);
                etaStripPeak.pk3  = getPeakBin20N(ecalRegion.etaStrip3);
                etaStripPeak.pk4  = getPeakBin20N(ecalRegion.etaStrip4);
                etaStripPeak.pk5  = getPeakBin20N(ecalRegion.etaStrip5);
                etaStripPeak.pk6  = getPeakBin20N(ecalRegion.etaStrip6);
                etaStripPeak.pk7  = getPeakBin20N(ecalRegion.etaStrip7);
                etaStripPeak.pk8  = getPeakBin20N(ecalRegion.etaStrip8);
                etaStripPeak.pk9  = getPeakBin20N(ecalRegion.etaStrip9);
                etaStripPeak.pk10 = getPeakBin20N(ecalRegion.etaStrip10);
                etaStripPeak.pk11 = getPeakBin20N(ecalRegion.etaStrip11);
                etaStripPeak.pk12 = getPeakBin20N(ecalRegion.etaStrip12);
                etaStripPeak.pk13 = getPeakBin20N(ecalRegion.etaStrip13);
                etaStripPeak.pk14 = getPeakBin20N(ecalRegion.etaStrip14);

        crystalMax peakIn15;
        peakIn15 = getPeakBin15N(etaStripPeak);

        cluster.seedEnergy 	= peakIn15.energy ;
        cluster.energy 		= 0 ;
        cluster.etaMax 		= peakIn15.etaMax ;
        cluster.phiMax 		= peakIn15.phiMax ;
        cluster.brems 		= 0 ;
        cluster.et5x5 		= 0 ;
        cluster.et2x5 		= 0 ;

return cluster ;
}

Cluster packCluster(ap_uint<15>& clusterEt, ap_uint<5>& etaMax_t, ap_uint<5>& phiMax_t, ap_uint<15>& et5x5_t, ap_uint<15>& et2x5_t, ap_uint<2>& brems_t ){

        ap_uint<12> peggedEt;
        Cluster pack;

        ap_uint<5>	towerEta 		= (etaMax_t)/5;
        ap_uint<2>	towerPhi 		= (phiMax_t)/5;
        ap_uint<3>	clusterEta 		= etaMax_t - 5*towerEta;
        ap_uint<3>	clusterPhi 		= phiMax_t - 5*towerPhi;
        ap_uint<15> clusterEt5x5 	= et5x5_t;
        ap_uint<15> clusterEt2x5 	= et2x5_t;
        ap_uint<2> 	clusterBrems 	= brems_t;

        peggedEt = (clusterEt > 0xFFF) ? (ap_uint<12>)0xFFF : (ap_uint<12>) clusterEt;

        pack = Cluster(peggedEt, towerEta, towerPhi, clusterEta, clusterPhi, 0, clusterEt5x5, clusterEt2x5, clusterBrems);

return pack;
}

void RemoveTmp(crystal temp[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> seed_eta,  ap_uint<5> seed_phi, ap_uint<2> brems  ){
#pragma HLS ARRAY_PARTITION variable=temp complete dim=0
#pragma HLS latency min=3
//#pragma HLS PIPELINE II=9

        for(loop i=0; i<CRYSTAL_IN_ETA; i++){
//                        #pragma HLS UNROLL
           if(i>=seed_eta-1 && i<=seed_eta+1){
           for(loop k=0; k<CRYSTAL_IN_PHI; k++){
                        #pragma HLS UNROLL
            if(k>=seed_phi-2 && k<=seed_phi+2)  temp[i][k].energy = 0 ;}
         }}

        if(brems == 1){
        for(loop i=0; i<CRYSTAL_IN_ETA; i++){
//                        #pragma HLS UNROLL
           if(i>=seed_eta-1 && i<=seed_eta+1){
           for(loop k=0; k<CRYSTAL_IN_PHI; k++){
                        #pragma HLS UNROLL
            if(k>=seed_phi-2-5 && k<=seed_phi+2-5)  temp[i][k].energy = 0 ;}
                        }
                }
         }

        if(brems == 2){
        for(loop i=0; i<CRYSTAL_IN_ETA; i++){
//                        #pragma HLS UNROLL
           if(i>=seed_eta-1 && i<=seed_eta+1 ){
           for(loop k=0; k<CRYSTAL_IN_PHI; k++){
                        #pragma HLS UNROLL
            if(k>=seed_phi-2+5 && k<=seed_phi+2+5)  temp[i][k].energy = 0 ;}
         }}}
}

clusterInfo getBremsValuesPos(crystal tempX[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> seed_eta,  ap_uint<5> seed_phi ){
#pragma HLS ARRAY_PARTITION variable=tempX complete dim=0
#pragma HLS latency min=6
//#pragma HLS PIPELINE II=9

        ap_uint<12> temp[CRYSTAL_IN_ETA+2][CRYSTAL_IN_PHI+4] ;
#pragma HLS ARRAY_PARTITION variable=temp complete dim=0

        ap_uint<12> eta_slice[3] ;
#pragma HLS ARRAY_PARTITION variable=eta_slice complete dim=0

clusterInfo cluster_tmp;

        for(loop i=0; i<CRYSTAL_IN_ETA+2; i++){
                        #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI+4; k++){
                        #pragma HLS UNROLL
            temp[i][k] = 0 ;
         }}

        for(loop i=0; i<CRYSTAL_IN_ETA; i++){
                        #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI-3; k++){
                        #pragma HLS UNROLL
            temp[i+1][k] = tempX[i][k+3].energy ;
         }}


        ap_uint<6> seed_eta1,  seed_phi1 ;

        seed_eta1 = seed_eta ; //to start from corner
        seed_phi1 = seed_phi ; //to start from corner
// now we are in the left bottom corner
        ap_uint<12> tmp1, tmp2, tmp3 ;

        for(loop j=0; j<CRYSTAL_IN_ETA; j++){
//                        #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI; k++){
                        #pragma HLS UNROLL
              if(j== seed_eta1 && k == seed_phi1)
                 {
                for(loop m=0; m<3 ; m++){
                        #pragma HLS UNROLL
                tmp1 = temp[j+m][k] + temp[j+m][k+1] ;
                tmp2 = temp[j+m][k+2] + temp[j+m][k+3] ;
                tmp3 = tmp1 + temp[j+m][k+4] ;
                eta_slice[m] = tmp2 + tmp3 ;
//                eta_slice[m] = temp[j+m][k] + temp[j+m][k+1] +temp[j+m][k+2] +temp[j+m][k+3] +temp[j+m][k+4] ;
                        }
               }
          }}

         cluster_tmp.energy=eta_slice[0] + eta_slice[1] + eta_slice[2] ;

return cluster_tmp ;
}


clusterInfo getBremsValuesNeg(crystal tempX[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> seed_eta,  ap_uint<5> seed_phi ){
#pragma HLS ARRAY_PARTITION variable=tempX complete dim=0
#pragma HLS latency min=6
//#pragma HLS PIPELINE II=9

        ap_uint<12> temp[CRYSTAL_IN_ETA+2][CRYSTAL_IN_PHI+4] ;
#pragma HLS ARRAY_PARTITION variable=temp complete dim=0

        ap_uint<12> eta_slice[3] ;
#pragma HLS ARRAY_PARTITION variable=eta_slice complete dim=0

clusterInfo cluster_tmp;

        for(loop i=0; i<CRYSTAL_IN_ETA+2; i++){
                        #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI+4; k++){
                        #pragma HLS UNROLL
            temp[i][k] = 0 ;
         }}

        for(loop i=0; i<CRYSTAL_IN_ETA; i++){
                        #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI-3; k++){
                        #pragma HLS UNROLL
            temp[i+1][k+7] = tempX[i][k].energy ;
         }}


        ap_uint<6> seed_eta1,  seed_phi1 ;

        seed_eta1 = seed_eta ; //to start from corner
        seed_phi1 = seed_phi ; //to start from corner
// now we are in the left bottom corner

        ap_uint<12> tmp1, tmp2, tmp3 ;

        for(loop j=0; j<CRYSTAL_IN_ETA; j++){
//                        #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI; k++){
                        #pragma HLS UNROLL
              if(j== seed_eta1 && k == seed_phi1)
                 {
                for(loop m=0; m<3 ; m++){
                        #pragma HLS UNROLL
                tmp1 = temp[j+m][k] + temp[j+m][k+1] ;
                tmp2 = temp[j+m][k+2] + temp[j+m][k+3] ;
                tmp3 = tmp1 + temp[j+m][k+4] ;
                eta_slice[m] = tmp2 + tmp3 ;
                        }
               }
          }}

         cluster_tmp.energy=eta_slice[0] + eta_slice[1] + eta_slice[2] ;

return cluster_tmp ;
}

clusterInfo getClusterValues(crystal tempX[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> seed_eta,  ap_uint<5> seed_phi ){
#pragma HLS ARRAY_PARTITION variable=tempX complete dim=0
#pragma HLS latency min=6
//#pragma HLS PIPELINE II=9

        ap_uint<12> temp[CRYSTAL_IN_ETA+4][CRYSTAL_IN_PHI+4] ;
#pragma HLS ARRAY_PARTITION variable=temp complete dim=0

        ap_uint<12> eta_slice[5] ;
#pragma HLS ARRAY_PARTITION variable=eta_slice complete dim=0


        ap_uint<12> et2x5_1Tot, et2x5_2Tot, etSum2x5 ;
        ap_uint<12> et5x5Tot ;

clusterInfo cluster_tmp;

        for(loop i=0; i<CRYSTAL_IN_ETA+4; i++){
                       #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI+4; k++){
                        #pragma HLS UNROLL
            temp[i][k] = 0 ;
         }}

        for(loop i=0; i<CRYSTAL_IN_ETA; i++){
                        #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI; k++){
                        #pragma HLS UNROLL
            temp[i+2][k+2] = tempX[i][k].energy ;
         }}


        ap_uint<6> seed_eta1,  seed_phi1 ;

        seed_eta1 = seed_eta ; //to start from corner
        seed_phi1 = seed_phi ; //to start from corner
// now we are in the left bottom corner
        ap_uint<12> tmp1, tmp2, tmp3 ;

        for(loop j=0; j<CRYSTAL_IN_ETA; j++){
//                        #pragma HLS UNROLL
           for(loop k=0; k<CRYSTAL_IN_PHI; k++){
                        #pragma HLS UNROLL
              if(j== seed_eta1 && k == seed_phi1)
                 {
                for(loop m=0; m<5 ; m++){
                        #pragma HLS UNROLL
                tmp1 = temp[j+m][k] + temp[j+m][k+1] ;
                tmp2 = temp[j+m][k+2] + temp[j+m][k+3] ;
                tmp3 = tmp1 + temp[j+m][k+4] ;
                eta_slice[m] = tmp2 + tmp3 ;
//                eta_slice[m] = temp[j+m][k] + temp[j+m][k+1] +temp[j+m][k+2] +temp[j+m][k+3] +temp[j+m][k+4] ;
                        }
               }
          }}


         cluster_tmp.energy=eta_slice[1] + eta_slice[2] + eta_slice[3] ;

          et5x5Tot = eta_slice[0] + eta_slice[1] + eta_slice[2] + eta_slice[3] + eta_slice[4] ;
          et2x5_1Tot = eta_slice[1] + eta_slice[2] ;
          et2x5_2Tot = eta_slice[2] + eta_slice[3] ;


          if(et2x5_1Tot >= et2x5_2Tot) etSum2x5 = et2x5_1Tot ;
          else etSum2x5 = et2x5_2Tot ;

          cluster_tmp.et5x5 = et5x5Tot ;
          cluster_tmp.et2x5 = etSum2x5 ;

return cluster_tmp ;
}


Cluster getRegion3x4(crystal temp[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> eta_offset){
#pragma HLS ARRAY_PARTITION variable=temp complete dim=0
//#pragma HLS PIPELINE II=2
#pragma HLS latency min=24

Cluster returnCluster, forCluster;
clusterInfo cluster_tmp;
clusterInfo cluster_tmpCenter;
clusterInfo cluster_tmpBneg;
clusterInfo cluster_tmpBpos;

ecalRegion_t ecalRegion;

        ecalRegion = initStructure(temp) ;

        cluster_tmp = getClusterPosition(ecalRegion) ;

        ap_uint<5> seed_phi = cluster_tmp.phiMax ;
        ap_uint<5> seed_eta = cluster_tmp.etaMax ;

        cluster_tmpCenter 	=  getClusterValues(temp, seed_eta, seed_phi) ;
        cluster_tmpBneg 	= getBremsValuesNeg(temp, seed_eta, seed_phi) ;
        cluster_tmpBpos 	= getBremsValuesPos(temp, seed_eta, seed_phi) ;

        cluster_tmp.energy = cluster_tmpCenter.energy;

        cluster_tmp.brems = 0 ;

        ap_uint<15> clusterEnergyDiv8 = 0;

        clusterEnergyDiv8 = cluster_tmpCenter.energy >> 3;

         if(cluster_tmpBneg.energy > clusterEnergyDiv8 && cluster_tmpBneg.energy > cluster_tmpBpos.energy) {
            cluster_tmp.energy = cluster_tmpCenter.energy + cluster_tmpBneg.energy;
            cluster_tmp.brems = 1;
         }
         else if(cluster_tmpBpos.energy > clusterEnergyDiv8){
            cluster_tmp.energy = cluster_tmpCenter.energy + cluster_tmpBpos.energy;
            cluster_tmp.brems  = 2;
         }

//eta, phi, seed energy in cluster_tmp; energy and brems in cluster_tmp1

        forCluster = packCluster(cluster_tmp.energy, cluster_tmp.etaMax, cluster_tmp.phiMax, cluster_tmp.et5x5, cluster_tmp.et2x5, cluster_tmp.brems);

        RemoveTmp(temp, seed_eta, seed_phi, cluster_tmp.brems);

        ap_uint<5> towerEta = forCluster.towerEta() + eta_offset;
        returnCluster = Cluster(forCluster.clusterEnergy(), towerEta, forCluster.towerPhi(), forCluster.clusterEta(), forCluster.clusterPhi(), forCluster.satur(), forCluster.et5x5(), forCluster.et2x5(), forCluster.brems());

return returnCluster;
}

void stitchClusters(Cluster ClusterUp, Cluster ClusterDown, Cluster &OutStchCluster1, Cluster &OutStchCluster2){
//#pragma HLS PIPELINE II=9

ap_uint<5> phi1 = ClusterUp.towerPhi()*5    + ClusterUp.clusterPhi();
ap_uint<5> phi2 = ClusterDown.towerPhi()*5  + ClusterDown.clusterPhi();
ap_uint<5> dPhi;

Cluster OutStchCluster1Reg = OutStchCluster1;
Cluster OutStchCluster2Reg = OutStchCluster2;

dPhi = (phi1>phi2)?(phi1-phi2):(phi2-phi1);

bool ClsInTwrEtaLevel1      = (ClusterDown.towerEta() == (ClusterUp.towerEta() + 1));
bool ClsInTwrEtaLevel2      = (ClusterDown.clusterEta() == 0 && ClusterUp.clusterEta() == 4);
bool ShrSamePhi1            = (dPhi<2);
//bool GrtThanZero			= (ClusterUp.clusterEnergy() > 0 && ClusterDown.clusterEnergy() > 0);

bool stitch                 = (ClsInTwrEtaLevel1 && ClsInTwrEtaLevel2 && ShrSamePhi1);

ap_uint<14> cEtSum          = ClusterUp.clusterEnergy() + ClusterDown.clusterEnergy();
ap_uint<12> pegged_cEtSum   = (cEtSum > 0xFFF) ? (ap_uint<12>)0xFFF : (ap_uint<12>) cEtSum;

ap_uint<17> et5x5           = ClusterUp.et5x5() + ClusterDown.et5x5();
ap_uint<15> pegged_et5x5    = (et5x5 > 0x7FFF) ? (ap_uint<15>)0x7FFF : (ap_uint<15>) et5x5;

ap_uint<17> et2x5           = ClusterUp.et2x5() + ClusterDown.et2x5();
ap_uint<15> pegged_et2x5    = (et2x5 > 0x7FFF) ? (ap_uint<15>)0x7FFF : (ap_uint<15>) et2x5;

if(stitch == 1){
    if(ClusterUp.clusterEnergy() > ClusterDown.clusterEnergy()){
                OutStchCluster1 = Cluster(pegged_cEtSum, ClusterUp.towerEta(), ClusterUp.towerPhi(), ClusterUp.clusterEta(), ClusterUp.clusterPhi(), ClusterUp.satur(), pegged_et5x5, pegged_et2x5, ClusterUp.brems());
                OutStchCluster2 = Cluster(0, 0, 0, 0, 0, 0, 0, 0, 0);
            }
            else{
                OutStchCluster1 = Cluster(0, 0, 0, 0, 0, 0, 0, 0, 0);
                OutStchCluster2 = Cluster(pegged_cEtSum, ClusterDown.towerEta(), ClusterDown.towerPhi(), ClusterDown.clusterEta(), ClusterDown.clusterPhi(), ClusterDown.satur(), pegged_et5x5, pegged_et2x5, ClusterDown.brems());
        }
}
else{
	OutStchCluster1 = OutStchCluster1Reg;
	OutStchCluster2 = OutStchCluster2Reg;
}
}

void algo_top(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]){
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS PIPELINE II=9
#pragma HLS INTERFACE ap_ctrl_hs port=return

//#pragma HLS latency min=100

crystal ECALRegion3x4_1[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI];
#pragma HLS ARRAY_PARTITION variable=ECALRegion3x4_1 complete dim=0

crystal ECALRegion3x4_2[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI];
#pragma HLS ARRAY_PARTITION variable=ECALRegion3x4_2 complete dim=0

crystal ECALRegion3x4_3[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI];
#pragma HLS ARRAY_PARTITION variable=ECALRegion3x4_3 complete dim=0

//creating 2 15x20 crystals temporary and temporary1

processInputLinks(link_in, ECALRegion3x4_1, ECALRegion3x4_2, ECALRegion3x4_3);

Cluster sort_clusterInStitched[Nbclusters];
Cluster sort_clusterIn[Nbclusters];
Cluster sort_clusterOut[Nbclusters];

tower_t towerEt[36];

#pragma HLS ARRAY_PARTITION variable=sort_clusterIn complete dim=0
#pragma HLS ARRAY_PARTITION variable=sort_clusterInStitched complete dim=0
#pragma HLS ARRAY_PARTITION variable=sort_clusterOut complete dim=0
#pragma HLS ARRAY_PARTITION variable=towerEt complete dim=0

// 9 is eta offset, since this is second part of RCT card
//
        sort_clusterIn[0]  = getRegion3x4(ECALRegion3x4_1, 9);
        sort_clusterIn[1]  = getRegion3x4(ECALRegion3x4_1, 9);
        sort_clusterIn[2]  = getRegion3x4(ECALRegion3x4_1, 9);
        sort_clusterIn[3]  = getRegion3x4(ECALRegion3x4_1, 9);

        sort_clusterIn[4]  = getRegion3x4(ECALRegion3x4_2, 12);
        sort_clusterIn[5]  = getRegion3x4(ECALRegion3x4_2, 12);
        sort_clusterIn[6]  = getRegion3x4(ECALRegion3x4_2, 12);
        sort_clusterIn[7]  = getRegion3x4(ECALRegion3x4_2, 12);

        sort_clusterIn[8]  = getRegion3x4(ECALRegion3x4_3, 15);
        sort_clusterIn[9]  = getRegion3x4(ECALRegion3x4_3, 15);
        sort_clusterIn[10] = getRegion3x4(ECALRegion3x4_3, 15);
        sort_clusterIn[11] = getRegion3x4(ECALRegion3x4_3, 15);

// stitch clusters

        for(loop i=0; i<12; i++){
                #pragma HLS UNROLL
                	sort_clusterInStitched[i] = sort_clusterIn[i];
                }

                for(loop i=12; i<16; i++){
                	#pragma HLS UNROLL
                	sort_clusterInStitched[i] = Cluster(0,0,0,0,0,0,0,0,0);
                }

                stitchClusters(sort_clusterIn[0], sort_clusterIn[4], sort_clusterInStitched[0], sort_clusterInStitched[4]);
                stitchClusters(sort_clusterIn[0], sort_clusterIn[5], sort_clusterInStitched[0], sort_clusterInStitched[5]);
                stitchClusters(sort_clusterIn[0], sort_clusterIn[6], sort_clusterInStitched[0], sort_clusterInStitched[6]);
                stitchClusters(sort_clusterIn[0], sort_clusterIn[7], sort_clusterInStitched[0], sort_clusterInStitched[7]);

                stitchClusters(sort_clusterIn[1], sort_clusterIn[4], sort_clusterInStitched[1], sort_clusterInStitched[4]);
                stitchClusters(sort_clusterIn[1], sort_clusterIn[5], sort_clusterInStitched[1], sort_clusterInStitched[5]);
                stitchClusters(sort_clusterIn[1], sort_clusterIn[6], sort_clusterInStitched[1], sort_clusterInStitched[6]);
                stitchClusters(sort_clusterIn[1], sort_clusterIn[7], sort_clusterInStitched[1], sort_clusterInStitched[7]);

                stitchClusters(sort_clusterIn[2], sort_clusterIn[4], sort_clusterInStitched[2], sort_clusterInStitched[4]);
                stitchClusters(sort_clusterIn[2], sort_clusterIn[5], sort_clusterInStitched[2], sort_clusterInStitched[5]);
                stitchClusters(sort_clusterIn[2], sort_clusterIn[6], sort_clusterInStitched[2], sort_clusterInStitched[6]);
                stitchClusters(sort_clusterIn[2], sort_clusterIn[7], sort_clusterInStitched[2], sort_clusterInStitched[7]);

                stitchClusters(sort_clusterIn[3], sort_clusterIn[4], sort_clusterInStitched[3], sort_clusterInStitched[4]);
                stitchClusters(sort_clusterIn[3], sort_clusterIn[5], sort_clusterInStitched[3], sort_clusterInStitched[5]);
                stitchClusters(sort_clusterIn[3], sort_clusterIn[6], sort_clusterInStitched[3], sort_clusterInStitched[6]);
                stitchClusters(sort_clusterIn[3], sort_clusterIn[7], sort_clusterInStitched[3], sort_clusterInStitched[7]);

                stitchClusters(sort_clusterIn[4], sort_clusterIn[8],  sort_clusterInStitched[4], sort_clusterInStitched[8]);
                stitchClusters(sort_clusterIn[4], sort_clusterIn[9],  sort_clusterInStitched[4], sort_clusterInStitched[9]);
                stitchClusters(sort_clusterIn[4], sort_clusterIn[10], sort_clusterInStitched[4], sort_clusterInStitched[10]);
                stitchClusters(sort_clusterIn[4], sort_clusterIn[11], sort_clusterInStitched[4], sort_clusterInStitched[11]);

                stitchClusters(sort_clusterIn[5], sort_clusterIn[8],  sort_clusterInStitched[5], sort_clusterInStitched[8]);
                stitchClusters(sort_clusterIn[5], sort_clusterIn[9],  sort_clusterInStitched[5], sort_clusterInStitched[9]);
                stitchClusters(sort_clusterIn[5], sort_clusterIn[10], sort_clusterInStitched[5], sort_clusterInStitched[10]);
                stitchClusters(sort_clusterIn[5], sort_clusterIn[11], sort_clusterInStitched[5], sort_clusterInStitched[11]);

                stitchClusters(sort_clusterIn[6], sort_clusterIn[8],  sort_clusterInStitched[6], sort_clusterInStitched[8]);
                stitchClusters(sort_clusterIn[6], sort_clusterIn[9],  sort_clusterInStitched[6], sort_clusterInStitched[9]);
                stitchClusters(sort_clusterIn[6], sort_clusterIn[10], sort_clusterInStitched[6], sort_clusterInStitched[10]);
                stitchClusters(sort_clusterIn[6], sort_clusterIn[11], sort_clusterInStitched[6], sort_clusterInStitched[11]);

                stitchClusters(sort_clusterIn[7], sort_clusterIn[8],  sort_clusterInStitched[7], sort_clusterInStitched[8]);
                stitchClusters(sort_clusterIn[7], sort_clusterIn[9],  sort_clusterInStitched[7], sort_clusterInStitched[9]);
                stitchClusters(sort_clusterIn[7], sort_clusterIn[10], sort_clusterInStitched[7], sort_clusterInStitched[10]);
                stitchClusters(sort_clusterIn[7], sort_clusterIn[11], sort_clusterInStitched[7], sort_clusterInStitched[11]);

// sorting of clusters, we have now 12 clusters, highest 6
// will be sent to SLRB, they come last from
// sorter

        bitonicSort16(sort_clusterInStitched, sort_clusterOut);

//for(loop i; i<16; i++) cout << sort_clusterOut[i].clusterEnergy() << endl ;
// build ECAL towers with unclustered energy
// keep only 6 highest clusters, others return back to towers

        ap_uint<12> towerEtECAL1[12];
        ap_uint<12> towerEtECAL2[12];
        ap_uint<12> towerEtECAL3[12];
        #pragma HLS ARRAY_PARTITION variable=towerEtECAL1 complete dim=0
        #pragma HLS ARRAY_PARTITION variable=towerEtECAL2 complete dim=0
        #pragma HLS ARRAY_PARTITION variable=towerEtECAL3 complete dim=0

        ap_uint<10> tower3x4_1_0[5][5];
                ap_uint<10> tower3x4_1_1[5][5];
                ap_uint<10> tower3x4_1_2[5][5];
                ap_uint<10> tower3x4_1_3[5][5];
                ap_uint<10> tower3x4_1_4[5][5];
                ap_uint<10> tower3x4_1_5[5][5];
                ap_uint<10> tower3x4_1_6[5][5];
                ap_uint<10> tower3x4_1_7[5][5];
                ap_uint<10> tower3x4_1_8[5][5];
                ap_uint<10> tower3x4_1_9[5][5];
                ap_uint<10> tower3x4_1_10[5][5];
                ap_uint<10> tower3x4_1_11[5][5];

                ap_uint<10> tower3x4_2_0[5][5];
                ap_uint<10> tower3x4_2_1[5][5];
                ap_uint<10> tower3x4_2_2[5][5];
                ap_uint<10> tower3x4_2_3[5][5];
                ap_uint<10> tower3x4_2_4[5][5];
                ap_uint<10> tower3x4_2_5[5][5];
                ap_uint<10> tower3x4_2_6[5][5];
                ap_uint<10> tower3x4_2_7[5][5];
                ap_uint<10> tower3x4_2_8[5][5];
                ap_uint<10> tower3x4_2_9[5][5];
                ap_uint<10> tower3x4_2_10[5][5];
                ap_uint<10> tower3x4_2_11[5][5];

                ap_uint<10> tower3x4_3_0[5][5];
                ap_uint<10> tower3x4_3_1[5][5];
                ap_uint<10> tower3x4_3_2[5][5];
                ap_uint<10> tower3x4_3_3[5][5];
                ap_uint<10> tower3x4_3_4[5][5];
                ap_uint<10> tower3x4_3_5[5][5];
                ap_uint<10> tower3x4_3_6[5][5];
                ap_uint<10> tower3x4_3_7[5][5];
                ap_uint<10> tower3x4_3_8[5][5];
                ap_uint<10> tower3x4_3_9[5][5];
                ap_uint<10> tower3x4_3_10[5][5];
                ap_uint<10> tower3x4_3_11[5][5];

                for(loop i=0; i<5; i++){
                	for(loop j=0; j<5; j++){
                		tower3x4_1_0[i][j]  = ECALRegion3x4_1[i+0][j+0].energy;
                		tower3x4_1_1[i][j]  = ECALRegion3x4_1[i+0][j+5].energy;
                		tower3x4_1_2[i][j]  = ECALRegion3x4_1[i+0][j+10].energy;
                		tower3x4_1_3[i][j]  = ECALRegion3x4_1[i+0][j+15].energy;
                		tower3x4_1_4[i][j]  = ECALRegion3x4_1[i+5][j+0].energy;
                		tower3x4_1_5[i][j]  = ECALRegion3x4_1[i+5][j+5].energy;
                		tower3x4_1_6[i][j]  = ECALRegion3x4_1[i+5][j+10].energy;
                		tower3x4_1_7[i][j]  = ECALRegion3x4_1[i+5][j+15].energy;
                		tower3x4_1_8[i][j]  = ECALRegion3x4_1[i+10][j+0].energy;
                		tower3x4_1_9[i][j]  = ECALRegion3x4_1[i+10][j+5].energy;
                		tower3x4_1_10[i][j] = ECALRegion3x4_1[i+10][j+10].energy;
                		tower3x4_1_11[i][j] = ECALRegion3x4_1[i+10][j+15].energy;

                		tower3x4_2_0[i][j]  = ECALRegion3x4_2[i+0][j+0].energy;
                		tower3x4_2_1[i][j]  = ECALRegion3x4_2[i+0][j+5].energy;
                		tower3x4_2_2[i][j]  = ECALRegion3x4_2[i+0][j+10].energy;
                		tower3x4_2_3[i][j]  = ECALRegion3x4_2[i+0][j+15].energy;
                		tower3x4_2_4[i][j]  = ECALRegion3x4_2[i+5][j+0].energy;
                		tower3x4_2_5[i][j]  = ECALRegion3x4_2[i+5][j+5].energy;
                		tower3x4_2_6[i][j]  = ECALRegion3x4_2[i+5][j+10].energy;
                		tower3x4_2_7[i][j]  = ECALRegion3x4_2[i+5][j+15].energy;
                		tower3x4_2_8[i][j]  = ECALRegion3x4_2[i+10][j+0].energy;
                		tower3x4_2_9[i][j]  = ECALRegion3x4_2[i+10][j+5].energy;
                		tower3x4_2_10[i][j] = ECALRegion3x4_2[i+10][j+10].energy;
                		tower3x4_2_11[i][j] = ECALRegion3x4_2[i+10][j+15].energy;

                		tower3x4_3_0[i][j]  = ECALRegion3x4_3[i+0][j+0].energy;
                		tower3x4_3_1[i][j]  = ECALRegion3x4_3[i+0][j+5].energy;
                		tower3x4_3_2[i][j]  = ECALRegion3x4_3[i+0][j+10].energy;
                		tower3x4_3_3[i][j]  = ECALRegion3x4_3[i+0][j+15].energy;
                		tower3x4_3_4[i][j]  = ECALRegion3x4_3[i+5][j+0].energy;
                		tower3x4_3_5[i][j]  = ECALRegion3x4_3[i+5][j+5].energy;
                		tower3x4_3_6[i][j]  = ECALRegion3x4_3[i+5][j+10].energy;
                		tower3x4_3_7[i][j]  = ECALRegion3x4_3[i+5][j+15].energy;
                		tower3x4_3_8[i][j]  = ECALRegion3x4_3[i+10][j+0].energy;
                		tower3x4_3_9[i][j]  = ECALRegion3x4_3[i+10][j+5].energy;
                		tower3x4_3_10[i][j] = ECALRegion3x4_3[i+10][j+10].energy;
                		tower3x4_3_11[i][j] = ECALRegion3x4_3[i+10][j+15].energy;
                	}
                }

                towerEtECAL1[0]  = getTowerEt(tower3x4_1_0);
                towerEtECAL1[1]  = getTowerEt(tower3x4_1_1);
                towerEtECAL1[2]  = getTowerEt(tower3x4_1_2);
                towerEtECAL1[3]  = getTowerEt(tower3x4_1_3);
                towerEtECAL1[4]  = getTowerEt(tower3x4_1_4);
                towerEtECAL1[5]  = getTowerEt(tower3x4_1_5);
                towerEtECAL1[6]  = getTowerEt(tower3x4_1_6);
                towerEtECAL1[7]  = getTowerEt(tower3x4_1_7);
                towerEtECAL1[8]  = getTowerEt(tower3x4_1_8);
                towerEtECAL1[9]  = getTowerEt(tower3x4_1_9);
                towerEtECAL1[10] = getTowerEt(tower3x4_1_10);
                towerEtECAL1[11] = getTowerEt(tower3x4_1_11);

                towerEtECAL2[0]  = getTowerEt(tower3x4_2_0);
                towerEtECAL2[1]  = getTowerEt(tower3x4_2_1);
                towerEtECAL2[2]  = getTowerEt(tower3x4_2_2);
                towerEtECAL2[3]  = getTowerEt(tower3x4_2_3);
                towerEtECAL2[4]  = getTowerEt(tower3x4_2_4);
                towerEtECAL2[5]  = getTowerEt(tower3x4_2_5);
                towerEtECAL2[6]  = getTowerEt(tower3x4_2_6);
                towerEtECAL2[7]  = getTowerEt(tower3x4_2_7);
                towerEtECAL2[8]  = getTowerEt(tower3x4_2_8);
                towerEtECAL2[9]  = getTowerEt(tower3x4_2_9);
                towerEtECAL2[10] = getTowerEt(tower3x4_2_10);
                towerEtECAL2[11] = getTowerEt(tower3x4_2_11);

                towerEtECAL3[0]  = getTowerEt(tower3x4_3_0);
                towerEtECAL3[1]  = getTowerEt(tower3x4_3_1);
                towerEtECAL3[2]  = getTowerEt(tower3x4_3_2);
                towerEtECAL3[3]  = getTowerEt(tower3x4_3_3);
                towerEtECAL3[4]  = getTowerEt(tower3x4_3_4);
                towerEtECAL3[5]  = getTowerEt(tower3x4_3_5);
                towerEtECAL3[6]  = getTowerEt(tower3x4_3_6);
                towerEtECAL3[7]  = getTowerEt(tower3x4_3_7);
                towerEtECAL3[8]  = getTowerEt(tower3x4_3_8);
                towerEtECAL3[9]  = getTowerEt(tower3x4_3_9);
                towerEtECAL3[10] = getTowerEt(tower3x4_3_10);
                towerEtECAL3[11] = getTowerEt(tower3x4_3_11);

        ap_uint<12> SumE1 ;
        ap_uint<12> SumE2 ;
        ap_uint<12> SumE3 ;

        ap_uint<12> SumE1c[6] ;
        ap_uint<12> SumE2c[6] ;
        ap_uint<12> SumE3c[6] ;
        #pragma HLS ARRAY_PARTITION variable=SumE1c complete dim=0
        #pragma HLS ARRAY_PARTITION variable=SumE2c complete dim=0
        #pragma HLS ARRAY_PARTITION variable=SumE3c complete dim=0

/* -9 takes into account eta offset of the clusters in this SLR */

        for(loop i=0; i<12; i++){
//1                #pragma HLS UNROLL
         for(loop k=4; k<10; k++){
                #pragma HLS UNROLL
          if(sort_clusterOut[k].clusterEnergy() > 0 && ((((sort_clusterOut[k].towerEta()-9)<<2)+sort_clusterOut[k].towerPhi()) == i)){
          SumE1c[k-4] = sort_clusterOut[k].clusterEnergy();
          sort_clusterOut[k] = Cluster(0, 0, 0, 0, 0, 0, 0, 0, 0);
          }
          else SumE1c[k-4] = 0;
          if(sort_clusterOut[k].clusterEnergy() > 0 && ((((sort_clusterOut[k].towerEta()-9)<<2)+sort_clusterOut[k].towerPhi()) == i+12)){
          SumE2c[k-4] = sort_clusterOut[k].clusterEnergy();
          sort_clusterOut[k] = Cluster(0, 0, 0, 0, 0, 0, 0, 0, 0);
          }
          else SumE2c[k-4] = 0;
          if(sort_clusterOut[k].clusterEnergy() > 0 && ((((sort_clusterOut[k].towerEta()-9)<<2)+sort_clusterOut[k].towerPhi()) == i+24)){
          SumE3c[k-4] = sort_clusterOut[k].clusterEnergy() ;
          sort_clusterOut[k] = Cluster(0, 0, 0, 0, 0, 0, 0, 0, 0);
          }
          else SumE3c[k-4] = 0;
         }
         	 SumE1 			= towerEtECAL1[i] + SumE1c[0] + SumE1c[1]+ SumE1c[2]+ SumE1c[3]+ SumE1c[4]+ SumE1c[5];
         	 SumE2 			= towerEtECAL2[i] + SumE2c[0] + SumE2c[1]+ SumE2c[2]+ SumE2c[3]+ SumE2c[4]+ SumE2c[5];
         	 SumE3 			= towerEtECAL3[i] + SumE3c[0] + SumE3c[1]+ SumE3c[2]+ SumE3c[3]+ SumE3c[4]+ SumE3c[5];
         	 towerEt[i]    	= tower_t(SumE1, 0, 0);
         	 towerEt[i+12] 	= tower_t(SumE2, 0, 0);
         	 towerEt[i+24] 	= tower_t(SumE3, 0, 0);
        }

        Cluster inClusterLink0[3];
                        Cluster inClusterLink1[3];

                        tower_t inTowerLink0[16];
                        tower_t inTowerLink1[16];

                        #pragma HLS ARRAY_PARTITION variable=inClusterLink0 complete dim=0
                        #pragma HLS ARRAY_PARTITION variable=inClusterLink0 complete dim=0
                        #pragma HLS ARRAY_PARTITION variable=inTowerLink0 complete dim=0
                        #pragma HLS ARRAY_PARTITION variable=inTowerLink1 complete dim=0



                        for(loop oLink=0; oLink<16; oLink++){
                                #pragma HLS unroll
                                inTowerLink0[oLink] = towerEt[oLink];
                                inTowerLink1[oLink] = towerEt[oLink+16];
                        }

                        for(loop oLink=0; oLink<3; oLink++){
                                #pragma HLS unroll
                                inClusterLink0[oLink] = sort_clusterOut[oLink+10];
                                inClusterLink1[oLink] = sort_clusterOut[oLink+13];
                        }

                        /*---------------------------------link 1------------------------------------*/
                         processOutLink(inClusterLink0, inTowerLink0, link_out[0]);

                        /*---------------------------------link 2------------------------------------*/
                         processOutLink(inClusterLink1, inTowerLink1, link_out[1]);

}

