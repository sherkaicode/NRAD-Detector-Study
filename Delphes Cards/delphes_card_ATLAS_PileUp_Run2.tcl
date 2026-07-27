#######################################
# Order of execution of various modules
#######################################
# Improved updated the parametric fit from the delphes_ATLAS_PileUp Card 
# Presentation, data, and code for the fits of the updates can be seen here: https://zenodo.org/records/17622130

set ExecutionPath {
    PileUpMerger
    ParticlePropagator

    ChargedHadronTrackingEfficiency
    ElectronTrackingEfficiency
    MuonTrackingEfficiency

    ChargedHadronMomentumSmearing
    ElectronMomentumSmearing
    MuonMomentumSmearing

    TrackMerger
    TrackSmearing
    Calorimeter
    ElectronFilter
    TrackPileUpSubtractor
    NeutralTowerMerger
    EFlowMergerAllTracks
    EFlowMerger
    EFlowFilter
    
    NeutrinoFilter
    GenJetFinder
    GenMissingET
    
    Rho
    FastJetFinder
    FatJetFinder
    JetPileUpSubtractor
    JetEnergyScale

    PhotonEfficiency
    PhotonIsolation

    ElectronEfficiency
    ElectronIsolation

    MuonEfficiency
    MuonIsolation

    MissingET
    
    JetFlavorAssociation

    BTagging
    TauTagging

    UniqueObjectFinder

    ScalarHT

    TreeWriter
}

#####################
# PileUp Merger
#####################

module PileUpMerger PileUpMerger {
    set InputArray Delphes/stableParticles

    set ParticleOutputArray stableParticles
    set VertexOutputArray vertices

    # pre-generated minbias input file
    set PileUpFile "/home/aegis/Titan0/madgraph/MG5_aMC_v3_7_0/bin/ttbar_validation/Events/Delphes Cards/13TeV_MinBias.pileup"

    # average expected pile up
    set MeanPileUp 40

    # maximum spread in the beam direction in m
    set ZVertexSpread 0.25

    # maximum spread in time in s
    set TVertexSpread 800E-12

    # vertex smearing formula f(z,t) (z,t need to be respectively given in m,s)
    set VertexDistributionFormula {exp(-(t^2/160e-12^2/2))*exp(-(z^2/0.053^2/2))}
}


#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {
    set InputArray PileUpMerger/stableParticles

    set OutputArray stableParticles
    set ChargedHadronOutputArray chargedHadrons
    set ElectronOutputArray electrons
    set MuonOutputArray muons

    set Radius 1.15
    set HalfLength 3.51

    set Bz 2.0
}

####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
    set InputArray ParticlePropagator/chargedHadrons
    set OutputArray chargedHadrons

    set EfficiencyFormula { 
        (abs(eta) <= 1.2) * ( 
            (pt < 200) * 0.95 + 
            ((pt >= 200) && (pt < 400)) * 0.93687 + 
            ((pt >= 400) && (pt < 600)) * 0.93296 + 
            ((pt >= 600) && (pt < 800)) * 0.92654 + 
            ((pt >= 800) && (pt < 1000)) * 0.92011 +
            ((pt >= 1000) && (pt < 1200)) * 0.91257 + 
            ((pt >= 1200) && (pt < 1400)) * 0.90531 + 
            ((pt >= 1400) && (pt < 1600)) * 0.89721 
        ) + 
        ((abs(eta) > 1.2) && (abs(eta) < 2.5)) * ( 
            (pt < 200) * 0.95 + 
            ((pt >= 200) && (pt < 400)) * 0.86145 + 
            ((pt >= 400) && (pt < 600)) * 0.85838 + 
            ((pt >= 600) && (pt < 800)) * 0.85251 + 
            ((pt >= 800) && (pt < 1000)) * 0.84749 + 
            ((pt >= 1000) && (pt < 1200)) * 0.83939 + 
            ((pt >= 1200) && (pt < 1400)) * 0.83212 + 
            ((pt >= 1400) && (pt < 1600)) * 0.82039 
        ) + 
        (abs(eta) >= 2.5) * 0.0 
    }

}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
    set InputArray ParticlePropagator/electrons
    set OutputArray electrons

    set EfficiencyFormula { 
        (pt <= 0.5) * (0.00) + 
        (abs(eta) <= 2.5) * ( 
            ((pt > 0.5) && (pt <= 10.0)) * 0.950 + 
            ((pt > 10.0) && (pt <= 15.0)) * 0.970 + 
            ((pt > 15.0) && (pt <= 30.0)) * 0.980 + 
            ((pt > 30.0) && (pt <= 50.0)) * 0.985 + 
            (pt > 50.0) * 0.990 
        ) + 
        (abs(eta) > 2.5) * 0.00 
    }
}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
    set InputArray ParticlePropagator/muons
    set OutputArray muons

    set EfficiencyFormula {
        (pt <= 0.1) * (0.00) +
        (abs(eta) <= 2.5) * (
            ((pt > 0.1) && (pt <= 5.0)) * 0.980 +
            (pt > 5.0) * 0.995
        ) +
        (abs(eta) > 2.5) * 0.00
    }
}

########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
    set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
    set OutputArray chargedHadrons
    set ResolutionFormula {                  (abs(eta) <= 0.5) * (pt > 0.1) * sqrt(0.06^2 + pt^2*1.3e-3^2) +
                            (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 0.1) * sqrt(0.10^2 + pt^2*1.7e-3^2) +
                            (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1) * sqrt(0.25^2 + pt^2*3.1e-3^2)
    }
}

###################################
# Momentum resolution for electrons
###################################

module MomentumSmearing ElectronMomentumSmearing {
    set InputArray ElectronTrackingEfficiency/electrons
    set OutputArray electrons
    set ResolutionFormula {
        (-2.47 < eta) && (eta <= -2.30) * (pt > 0.1) * sqrt((0.025/sqrt(pt))^2 + (0.02544753)^2) +
        (-2.30 < eta) && (eta <= -2.00) * (pt > 0.1) * sqrt((0.009/sqrt(pt))^2 + (0.01294753)^2) +
        (-2.00 < eta) && (eta <= -1.80) * (pt > 0.1) * sqrt((0.006/sqrt(pt))^2 + (0.01271605)^2) +
        (-1.80 < eta) && (eta <= -1.60) * (pt > 0.1) * sqrt((0.005/sqrt(pt))^2 + (0.01557099)^2) +
        (-1.60 < eta) && (eta <= -1.40) * (pt > 0.1) * sqrt((0.011/sqrt(pt))^2 + (0.02444444)^2) +
        (-1.40 < eta) && (eta <= -1.20) * (pt > 0.1) * sqrt((0.002/sqrt(pt))^2 + (0.01595679)^2) +
        (-1.20 < eta) && (eta <= -1.00) * (pt > 0.1) * sqrt((0.017/sqrt(pt))^2 + (0.00677469)^2) +
        (-1.00 < eta) && (eta <= -0.80) * (pt > 0.1) * sqrt((0.011/sqrt(pt))^2 + (0.00808642)^2) +
        (-0.80 < eta) && (eta <= -0.60) * (pt > 0.1) * sqrt((0.020/sqrt(pt))^2 + (0.00916667)^2) +
        (-0.60 < eta) && (eta <= -0.40) * (pt > 0.1) * sqrt((0.019/sqrt(pt))^2 + (0.00754630)^2) +
        (-0.40 < eta) && (eta <= -0.20) * (pt > 0.1) * sqrt((0.021/sqrt(pt))^2 + (0.00631173)^2) +
        (-0.20 < eta) && (eta <= 0.00) * (pt > 0.1) * sqrt((0.017/sqrt(pt))^2 + (0.00816358)^2) +
        (0.00 < eta) && (eta <= 0.20) * (pt > 0.1) * sqrt((0.015/sqrt(pt))^2 + (0.00754630)^2) +
        (0.20 < eta) && (eta <= 0.40) * (pt > 0.1) * sqrt((0.019/sqrt(pt))^2 + (0.00685185)^2) +
        (0.40 < eta) && (eta <= 0.60) * (pt > 0.1) * sqrt((0.018/sqrt(pt))^2 + (0.00777778)^2) +
        (0.60 < eta) && (eta <= 0.80) * (pt > 0.1) * sqrt((0.019/sqrt(pt))^2 + (0.00901235)^2) +
        (0.80 < eta) && (eta <= 1.00) * (pt > 0.1) * sqrt((0.011/sqrt(pt))^2 + (0.01202160)^2) +
        (1.00 < eta) && (eta <= 1.20) * (pt > 0.1) * sqrt((0.016/sqrt(pt))^2 + (0.01086420)^2) +
        (1.20 < eta) && (eta <= 1.40) * (pt > 0.1) * sqrt((0.002/sqrt(pt))^2 + (0.01564815)^2) +
        (1.40 < eta) && (eta <= 1.60) * (pt > 0.1) * sqrt((0.018/sqrt(pt))^2 + (0.01966049)^2) +
        (1.60 < eta) && (eta <= 1.80) * (pt > 0.1) * sqrt((0.001/sqrt(pt))^2 + (0.01580247)^2) +
        (1.80 < eta) && (eta <= 2.00) * (pt > 0.1) * sqrt((0.013/sqrt(pt))^2 + (0.01148148)^2) +
        (2.00 < eta) && (eta <= 2.30) * (pt > 0.1) * sqrt((0.004/sqrt(pt))^2 + (0.01433642)^2) +
        (2.30 < eta) && (eta <= 2.47) * (pt > 0.1) * sqrt((0.017/sqrt(pt))^2 + (0.02899691)^2)
    }
}

###############################
# Momentum resolution for muons 
###############################

module MomentumSmearing MuonMomentumSmearing {
    set InputArray MuonTrackingEfficiency/muons
    set OutputArray muons
    set ResolutionFormula { 
        (abs(eta) <= 1.05) * sqrt((((6.7+6.5)/2.0*0.001)^2) + (((0.08+0.11)/2.0*0.001*pt)^2)) + 
        (abs(eta) > 1.05 && abs(eta) < 2.0) * sqrt((((10.3+8.9)/2.0*0.001)^2) + (((0.24+0.29)/2.0*0.001*pt)^2)) +
        (abs(eta) >= 2.0) * sqrt((((10.6+11.5)/2.0*0.001)^2) +(((0.21+0.26)/2.0*0.001*pt)^2))
    }
}

##############
# Track merger
##############
  
module Merger TrackMerger {
    add InputArray ChargedHadronMomentumSmearing/chargedHadrons
    add InputArray ElectronMomentumSmearing/electrons
    add InputArray MuonMomentumSmearing/muons
    set OutputArray tracks
}
module TrackSmearing TrackSmearing {
    set InputArray TrackMerger/tracks
    set OutputArray tracks
    set Bz 2.0

    source trackResolutionATLAS.tcl
}

#############
# Calorimeter
#############

module Calorimeter Calorimeter {
    set ParticleInputArray ParticlePropagator/stableParticles
    set TrackInputArray TrackSmearing/tracks
    set TowerOutputArray towers
    set PhotonOutputArray photons

    set EFlowTrackOutputArray eflowTracks
    set EFlowPhotonOutputArray eflowPhotons
    set EFlowNeutralHadronOutputArray eflowNeutralHadrons

    set ECalEnergyMin 0.5
    set HCalEnergyMin 1.0

    set ECalEnergySignificanceMin 1.0
    set HCalEnergySignificanceMin 1.0

    set SmearTowerCenter true

    set pi [expr {acos(-1)}]

    # 10 degrees towers
    set PhiBins {}
    for {set i -18} {$i <= 18} {incr i} {
        add PhiBins [expr {$i * $pi/18.0}]
    }
    foreach eta {-3.2 -2.5 -2.4 -2.3 -2.2 -2.1 -2 -1.9 -1.8 -1.7 -1.6 -1.5 -1.4 -1.3 -1.2 -1.1 -1 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 3.3} {
        add EtaPhiBins $eta $PhiBins
    }

    # 20 degrees towers
    set PhiBins {}
    for {set i -9} {$i <= 9} {incr i} {
        add PhiBins [expr {$i * $pi/9.0}]
    }
    foreach eta {-4.9 -4.7 -4.5 -4.3 -4.1 -3.9 -3.7 -3.5 -3.3 -3 -2.8 -2.6 2.8 3 3.2 3.5 3.7 3.9 4.1 4.3 4.5 4.7 4.9} {
        add EtaPhiBins $eta $PhiBins
    }

    # default energy fractions {abs(PDG code)} {Fecal Fhcal}
    add EnergyFraction {0} {0.0 1.0}
    # energy fractions for e, gamma and pi0
    add EnergyFraction {11} {1.0 0.0}
    add EnergyFraction {22} {1.0 0.0}
    add EnergyFraction {111} {1.0 0.0}
    # energy fractions for muon, neutrinos and neutralinos
    add EnergyFraction {12} {0.0 0.0}
    add EnergyFraction {13} {0.0 0.0}
    add EnergyFraction {14} {0.0 0.0}
    add EnergyFraction {16} {0.0 0.0}
    add EnergyFraction {1000022} {0.0 0.0}
    add EnergyFraction {1000023} {0.0 0.0}
    add EnergyFraction {1000025} {0.0 0.0}
    add EnergyFraction {1000035} {0.0 0.0}
    add EnergyFraction {1000045} {0.0 0.0}

    # --- BSM HV INVISIBLES ---
    # Visible Diagonals (Included here in case they are set to be stable)
    add EnergyFraction {4900111} {0.0 0.0}
    add EnergyFraction {4900113} {0.0 0.0}

    # Simple Setup Stable Hadrons (Off-Diagonals & Glueballs)
    add EnergyFraction {4900211} {0.0 0.0}
    add EnergyFraction {-4900211} {0.0 0.0}
    add EnergyFraction {4900213} {0.0 0.0}
    add EnergyFraction {-4900213} {0.0 0.0}
    add EnergyFraction {4900991} {0.0 0.0}

    # Hidden Baryons (Deltav)
    add EnergyFraction {4901114} {0.0 0.0}
    add EnergyFraction {-4901114} {0.0 0.0}

    # Hidden Quarks (Stable in U(1) or pre-hadronization)
    add EnergyFraction {4900101} {0.0 0.0}
    add EnergyFraction {4900102} {0.0 0.0}
    add EnergyFraction {4900103} {0.0 0.0}
    add EnergyFraction {4900104} {0.0 0.0}
    add EnergyFraction {4900105} {0.0 0.0}
    add EnergyFraction {4900106} {0.0 0.0}
    add EnergyFraction {4900107} {0.0 0.0}
    add EnergyFraction {4900108} {0.0 0.0}

    # Hidden Gauge Bosons
    add EnergyFraction {4900021} {0.0 0.0} 
    add EnergyFraction {4900022} {0.0 0.0} 

    # Generic Dark Matter (The ones used in your specific Rinv decay)
    add EnergyFraction {51} {0.0 0.0}
    add EnergyFraction {-51} {0.0 0.0}
    add EnergyFraction {53} {0.0 0.0}
    add EnergyFraction {-53} {0.0 0.0}
    add EnergyFraction {52} {0.0 0.0}
    add EnergyFraction {54} {0.0 0.0}

    add EnergyFraction {4900121} {0.0 0.0}
    add EnergyFraction {4900123} {0.0 0.0}
    add EnergyFraction {4900231} {0.0 0.0}
    add EnergyFraction {4900233} {0.0 0.0}
    
    # energy fractions for K0short and Lambda
    add EnergyFraction {310} {0.3 0.7}
    add EnergyFraction {3122} {0.3 0.7}

    # set ECalResolutionFormula {resolution formula as a function of eta and energy}
    # http://arxiv.org/pdf/physics/0608012v1 jinst8_08_s08003
    # http://villaolmo.mib.infn.it/ICATPP9th_2005/Calorimetry/Schram.p.pdf
    # http://www.physics.utoronto.ca/~krieger/procs/ComoProceedings.pdf
    set ECalResolutionFormula {                  (abs(eta) <= 3.2) * sqrt(energy^2*0.0017^2 + energy*0.101^2) +
                                (abs(eta) > 3.2 && abs(eta) <= 4.9) * sqrt(energy^2*0.0350^2 + energy*0.285^2)
    }

    # set HCalResolutionFormula {resolution formula as a function of eta and energy}
    # http://arxiv.org/pdf/hep-ex/0004009v1
    # http://villaolmo.mib.infn.it/ICATPP9th_2005/Calorimetry/Schram.p.pdf
    set HCalResolutionFormula {                  (abs(eta) <= 1.7) * sqrt(energy^2*0.0302^2 + energy*0.5205^2 + 1.59^2) +
                                (abs(eta) > 1.7 && abs(eta) <= 3.2) * sqrt(energy^2*0.0500^2 + energy*0.706^2) +
                                (abs(eta) > 3.2 && abs(eta) <= 4.9) * sqrt(energy^2*0.09420^2 + energy*1.00^2)                            
    }
}

#################
# Electron filter
#################

module PdgCodeFilter ElectronFilter {
    set InputArray Calorimeter/eflowTracks
    set OutputArray electrons
    set Invert true
    add PdgCode {11}
    add PdgCode {-11}
}

##########################
# Track pile-up subtractor
##########################

module TrackPileUpSubtractor TrackPileUpSubtractor {
    add InputArray Calorimeter/eflowTracks eflowTracks
    add InputArray ElectronFilter/electrons electrons
    add InputArray MuonMomentumSmearing/muons muons

    set VertexInputArray PileUpMerger/vertices
    set ZVertexResolution {0.0001}
}

####################
# Neutral tower merger
####################

module Merger NeutralTowerMerger {
    add InputArray Calorimeter/eflowPhotons
    add InputArray Calorimeter/eflowNeutralHadrons
    set OutputArray eflowTowers
}

##################################
# Energy flow merger (all tracks)
##################################

module Merger EFlowMergerAllTracks {
    add InputArray TrackSmearing/tracks
    add InputArray Calorimeter/eflowPhotons
    add InputArray Calorimeter/eflowNeutralHadrons
    set OutputArray eflow
}

####################
# Energy flow merger
####################

module Merger EFlowMerger {
    add InputArray Calorimeter/eflowTracks
    add InputArray Calorimeter/eflowPhotons
    add InputArray Calorimeter/eflowNeutralHadrons
    set OutputArray eflow
}

######################
# EFlowFilter
######################

module PdgCodeFilter EFlowFilter {
    set InputArray EFlowMerger/eflow
    set OutputArray eflow
    
    add PdgCode {11}
    add PdgCode {-11}
    add PdgCode {13}
    add PdgCode {-13}
}

#############
# Rho pile-up
#############

module FastJetGridMedianEstimator Rho {

    set InputArray Calorimeter/towers
    set RhoOutputArray rho

    add GridRange -5.0 -2.5 1.0 1.0
    add GridRange -2.5 2.5 0.5 0.5
    add GridRange 2.5 5.0 1.0 1.0

}

#####################
# Neutrino Filter
#####################

module PdgCodeFilter NeutrinoFilter {
    set InputArray Delphes/stableParticles
    set OutputArray filteredParticles

    set PTMin 0.0

    add PdgCode {12}
    add PdgCode {14}
    add PdgCode {16}
    add PdgCode {-12}
    add PdgCode {-14}
    add PdgCode {-16}
}

#####################
# MC truth jet finder
#####################

module FastJetFinder GenJetFinder {
    set InputArray NeutrinoFilter/filteredParticles

    set OutputArray jets

    # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
    set JetAlgorithm 6
    set ParameterR 0.4

    set JetPTMin 20.0
}

#########################
# Gen Missing ET merger
########################

module Merger GenMissingET {
    add InputArray NeutrinoFilter/filteredParticles
    set MomentumOutputArray momentum
}


############
# Jet finder
############

module FastJetFinder FastJetFinder {
    set InputArray Calorimeter/towers
    set OutputArray jets

    set AreaAlgorithm 5
    set JetAlgorithm 6
    set ParameterR 0.4
    set JetPTMin 20.0
}

##################
# Fat Jet finder
##################

module FastJetFinder FatJetFinder {
    set InputArray Calorimeter/towers
    set OutputArray jets
    set JetAlgorithm 6
    set ParameterR 1.0

    set ComputeNsubjettiness 1
    set Beta 1.0
    set AxisMode 4

    set ComputeTrimming 1
    set RTrim 0.2
    set PtFracTrim 0.05

    set JetPTMin 250.0
}

###########################
# Jet Pile-Up Subtraction
###########################

module JetPileUpSubtractor JetPileUpSubtractor {
    set JetInputArray FastJetFinder/jets
    set RhoInputArray Rho/rho

    set OutputArray jets

    set JetPTMin 20.0
}

##################
# Jet Energy Scale
##################

module EnergyScale JetEnergyScale {
    set InputArray JetPileUpSubtractor/jets
    set OutputArray jets
    set ScaleFormula {  sqrt( (3.5 - 0.2*(abs(eta)))^2 / pt + 1.0 )  }
}

###################
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
    set InputArray Calorimeter/eflowPhotons
    set OutputArray photons

    set EfficiencyFormula { 
        ((abs(eta) >= 0.0) && (abs(eta) < 0.6)) * ( 
            (pt > 25 && pt <= 30)  * 0.8379 + 
            (pt > 30 && pt <= 35)  * 0.8598 + 
            (pt > 35 && pt <= 40)  * 0.8766 + 
            (pt > 40 && pt <= 45)  * 0.8923 + 
            (pt > 45 && pt <= 50)  * 0.9038 + 
            (pt > 50 && pt <= 60)  * 0.9142 + 
            (pt > 60 && pt <= 80)  * 0.9257 + 
            (pt > 80 && pt <= 100) * 0.9320 + 
            (pt > 100 && pt <= 125) * 0.9362 + 
            (pt > 125 && pt <= 150) * 0.9383 + 
            (pt > 150 && pt <= 175) * 0.9383 + 
            (pt > 175 && pt <= 250) * 0.9351 + 
            (pt > 250 && pt <= 350) * 0.9331 + 
            (pt > 350 && pt <= 1000) * 0.9215 + 
            (pt > 1000) * 0.9215 
        ) + 
        ((abs(eta) >= 0.6) && (abs(eta) < 1.37)) * ( 
            (pt > 25 && pt <= 30)  * 0.8305 + 
            (pt > 30 && pt <= 35)  * 0.8515 + 
            (pt > 35 && pt <= 40)  * 0.8881 + 
            (pt > 40 && pt <= 45)  * 0.8828 + 
            (pt > 45 && pt <= 50)  * 0.8964 + 
            (pt > 50 && pt <= 60)  * 0.9090 + 
            (pt > 60 && pt <= 80)  * 0.9195 + 
            (pt > 80 && pt <= 100) * 0.9268 + 
            (pt > 100 && pt <= 125) * 0.9351 + 
            (pt > 125 && pt <= 150) * 0.9372 + 
            (pt > 150 && pt <= 175) * 0.9362 + 
            (pt > 175 && pt <= 250) * 0.9362 + 
            (pt > 250 && pt <= 350) * 0.9341 + 
            (pt > 350 && pt <= 1000) * 0.9299 + 
            (pt > 1000) * 0.9299 
        ) + 
        ((abs(eta) >= 1.52) && (abs(eta) < 1.81)) * ( 
            (pt > 25 && pt <= 30)  * 0.8421 + 
            (pt > 30 && pt <= 35)  * 0.8630 + 
            (pt > 35 && pt <= 40)  * 0.8577 + 
            (pt > 40 && pt <= 45)  * 0.8985 + 
            (pt > 45 && pt <= 50)  * 0.9100 + 
            (pt > 50 && pt <= 60)  * 0.9278 + 
            (pt > 60 && pt <= 80)  * 0.9351 + 
            (pt > 80 && pt <= 100) * 0.9393 + 
            (pt > 100 && pt <= 125) * 0.9404 + 
            (pt > 125 && pt <= 150) * 0.9414 + 
            (pt > 150 && pt <= 175) * 0.9446 + 
            (pt > 175 && pt <= 250) * 0.9467 + 
            (pt > 250 && pt <= 350) * 0.9446 + 
            (pt > 350 && pt <= 1000) * 0.9393 + 
            (pt > 1000) * 0.9393 
        ) + 
        ((abs(eta) >= 1.81) && (abs(eta) < 2.37)) * ( 
            (pt > 25 && pt <= 30)  * 0.8358 + 
            (pt > 30 && pt <= 35)  * 0.8483 + 
            (pt > 35 && pt <= 40)  * 0.8818 + 
            (pt > 40 && pt <= 45)  * 0.8808 + 
            (pt > 45 && pt <= 50)  * 0.8860 + 
            (pt > 50 && pt <= 60)  * 0.8923 + 
            (pt > 60 && pt <= 80)  * 0.8954 + 
            (pt > 80 && pt <= 100) * 0.9059 + 
            (pt > 100 && pt <= 125) * 0.9132 + 
            (pt > 125 && pt <= 150) * 0.9153 + 
            (pt > 150 && pt <= 175) * 0.9153 + 
            (pt > 175 && pt <= 250) * 0.9121 + 
            (pt > 250 && pt <= 350) * 0.9027 + 
            (pt > 350 && pt <= 1000) * 0.8923 + 
            (pt > 1000) * 0.8923 
        ) 
    }

}

##################
# Photon isolation
##################

module Isolation PhotonIsolation {
    set CandidateInputArray PhotonEfficiency/photons
    set IsolationInputArray EFlowFilter/eflow
    set RhoInputArray Rho/rho

    set OutputArray photons
    set DeltaRMax 0.2

    set PTMin 1.0
    set PTRatioMax 0.05
}

#####################
# Electron efficiency
#####################

module Efficiency ElectronEfficiency {
    add InputArray TrackPileUpSubtractor/electrons
    set OutputArray electrons
    set EfficiencyFormula { 
        (pt < 4.0 && abs(eta) < 2.47) * ( 
        (eta >= -2.47 && eta <= -2.35) * 0.7582 + 
        (eta > -2.35 && eta <= -2.00) * 0.8092 + 
        (eta > -2.00 && eta <= -1.80) * 0.8314 + 
        (eta > -1.80 && eta <= -1.50) * 0.8614 + 
        (eta > -1.50 && eta <= -1.35) * 0.8065 + 
        (eta > -1.35 && eta <= -1.15) * 0.8536 + 
        (eta > -1.15 && eta <= -0.80) * 0.8824 + 
        (eta > -0.80 && eta <= -0.60) * 0.8928 + 
        (eta > -0.60 && eta <= -0.01) * 0.8967 + 
        (eta > -0.01 && eta <= 0.00) * 0.8627 + 
        (eta > 0.00 && eta <= 0.01) * 0.8458 + 
        (eta > 0.01 && eta <= 0.60) * 0.8980 + 
        (eta > 0.60 && eta <= 0.80) * 0.8967 + 
        (eta > 0.80 && eta <= 1.15) * 0.8863 + 
        (eta > 1.15 && eta <= 1.35) * 0.8575 + 
        (eta > 1.35 && eta <= 1.50) * 0.8235 + 
        (eta > 1.50 && eta <= 1.80) * 0.8693 + 
        (eta > 1.80 && eta <= 2.00) * 0.8340 + 
        (eta > 2.00 && eta <= 2.35) * 0.8118 + 
        (eta > 2.35 && eta <= 2.47) * 0.7699) * ( 
        (pt > 4.0 && pt <= 7.3) * 0.7550 + 
        (pt > 7.3 && pt <= 10.0) * 0.7876 + 
        (pt > 10.0 && pt <= 15.0) * 0.7993 + 
        (pt > 15.0 && pt <= 20.0) * 0.7759 + 
        (pt > 20.0 && pt <= 25.0) * 0.7889 + 
        (pt > 25.0 && pt <= 30.0) * 0.8241 + 
        (pt > 30.0 && pt <= 35.0) * 0.8489 + 
        (pt > 35.0 && pt <= 40.0) * 0.8684 + 
        (pt > 40.0 && pt <= 45.0) * 0.8801 + 
        (pt > 45.0 && pt <= 50.0) * 0.8866 + 
        (pt > 50.0 && pt <= 60.0) * 0.8984 + 
        (pt > 60.0 && pt <= 80.0) * 0.9179 + 
        (pt > 80.0) * 0.9309) 
    }
}

####################
# Electron isolation
####################

module Isolation ElectronIsolation {
    set CandidateInputArray ElectronEfficiency/electrons
    set IsolationInputArray EFlowFilter/eflow
    set RhoInputArray Rho/rho

    set OutputArray electrons
    set DeltaRMax 0.2

    set PTMin 1.0
    set PTRatioMax 0.15
}

#################
# Muon efficiency
#################

module Efficiency MuonEfficiency {
    add InputArray TrackPileUpSubtractor/muons
    set OutputArray muons
    set EfficiencyFormula { 
        (pt > 3.0) * (abs(eta) < 2.5) * ( 
            ( 
                ((eta > -2.5) && (eta <= -2.279)) * 0.987 + 
                ((eta > -2.279) && (eta <= -2.047)) * 0.987 + 
                ((eta > -2.047) && (eta <= -1.82)) * 0.984 + 
                ((eta > -1.82) && (eta <= -1.589)) * 0.99 + 
                ((eta > -1.589) && (eta <= -1.367)) * 0.991 + 
                ((eta > -1.367) && (eta <= -1.135)) * 0.98 + 
                ((eta > -1.135) && (eta <= -0.909)) * 0.979 + 
                ((eta > -0.909) && (eta <= -0.682)) * 0.987 + 
                ((eta > -0.682) && (eta <= -0.456)) * 0.99 + 
                ((eta > -0.456) && (eta <= -0.229)) * 0.987 + 
                ((eta > -0.229) && (eta <= 0.003)) * 0.842 + 
                ((eta > 0.003) && (eta <= 0.229)) * 0.852 + 
                ((eta > 0.229) && (eta <= 0.456)) * 0.986 + 
                ((eta > 0.456) && (eta <= 0.682)) * 0.99 + 
                ((eta > 0.682) && (eta <= 0.909)) * 0.99 + 
                ((eta > 0.909) && (eta <= 1.141)) * 0.98 + 
                ((eta > 1.141) && (eta <= 1.362)) * 0.98 + 
                ((eta > 1.362) && (eta <= 1.594)) * 0.992 + 
                ((eta > 1.594) && (eta <= 1.82)) * 0.99 + 
                ((eta > 1.82) && (eta <= 2.047)) * 0.984 + 
                ((eta > 2.047) && (eta <= 2.273)) * 0.988 + 
                ((eta > 2.273) && (eta <= 2.5)) * 0.988 
            ) * 
            ( 
                ((pt > 3.015) && (pt <= 3.499)) * 0.453 + 
                ((pt > 3.499) && (pt <= 4.0)) * 0.614 + 
                ((pt > 4.0) && (pt <= 5.003)) * 0.808 + 
                ((pt > 5.003) && (pt <= 6.005)) * 0.94 + 
                ((pt > 6.005) && (pt <= 7.008)) * 0.963 + 
                ((pt > 7.008) && (pt <= 8.01)) * 0.973 + 
                ((pt > 8.01) && (pt <= 8.995)) * 0.972 + 
                ((pt > 8.995) && (pt <= 10.015)) * 0.975 + 
                ((pt > 10.015) && (pt <= 12.003)) * 0.975 + 
                ((pt > 12.003) && (pt <= 14.008)) * 0.975 + 
                ((pt > 14.008) && (pt <= 16.013)) * 0.975 + 
                ((pt > 16.013) && (pt <= 18.0)) * 0.975 + 
                ((pt > 18.0) && (pt <= 20.005)) * 0.977 + 
                (pt >= 20.005) * 0.977 
            ) 
        ) 
    }

}
################
# Muon isolation
################

module Isolation MuonIsolation {
    set CandidateInputArray MuonEfficiency/muons
    set IsolationInputArray EFlowFilter/eflow
    set RhoInputArray Rho/rho

    set OutputArray muons
    set DeltaRMax 0.3

    set PTMin 1.0
    set PTRatioMax 0.3
}

###################
# Missing ET merger
###################

module Merger MissingET {
    add InputArray EFlowMergerAllTracks/eflow
    set MomentumOutputArray momentum
}


##################
# Scalar HT merger
##################

module Merger ScalarHT {
    add InputArray UniqueObjectFinder/jets
    add InputArray UniqueObjectFinder/electrons
    add InputArray UniqueObjectFinder/photons
    add InputArray UniqueObjectFinder/muons
    set EnergyOutputArray energy
}

########################
# Jet Flavor Association
########################

module JetFlavorAssociation JetFlavorAssociation {
    set PartonInputArray Delphes/partons
    set ParticleInputArray Delphes/allParticles
    set ParticleLHEFInputArray Delphes/allParticlesLHEF
    set JetInputArray JetEnergyScale/jets
    
    set DeltaR 0.5
    set PartonPTMin 1.0
    set PartonEtaMax 2.5
}

###########
# b-tagging
###########

module BTagging BTagging {
    set JetInputArray JetEnergyScale/jets

    set BitNumber 0

    add EfficiencyFormula {0} {0.00677 + 2.1e-06*pt}
    add EfficiencyFormula {4} {0.186*tanh(0.60700*pt)*(1/(1 + 0.00097*pt))}
    add EfficiencyFormula {5} {2.993*tanh(0.00181*pt)*(30/(1 + 0.18066*pt))}
}

#############
# tau-tagging
#############

module TrackCountingTauTagging TauTagging {
    set ParticleInputArray Delphes/allParticles
    set PartonInputArray Delphes/partons
    set TrackInputArray TrackSmearing/tracks
    set JetInputArray JetEnergyScale/jets

    set DeltaR 0.2
    set DeltaRTrack 0.2

    set TrackPTMin 1.0
    
    set TauPTMin 1.0
    set TauEtaMax 2.5
    
    set BitNumber 0
    
    add EfficiencyFormula {1} {0.70}
    add EfficiencyFormula {2} {0.60}
    add EfficiencyFormula {-1} {0.02}
    add EfficiencyFormula {-2} {0.01}

}

#####################################################
# Find uniquely identified photons/electrons/tau/jets
#####################################################

module UniqueObjectFinder UniqueObjectFinder {
    add InputArray PhotonIsolation/photons photons
    add InputArray ElectronIsolation/electrons electrons
    add InputArray MuonIsolation/muons muons
    add InputArray JetEnergyScale/jets jets
}

##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
    add Branch Delphes/allParticles Particle GenParticle
    add Branch Calorimeter/towers Tower Tower
    add Branch EFlowMerger/eflow ParticleFlowCandidate ParticleFlowCandidate
    add Branch GenJetFinder/jets GenJet Jet
    add Branch GenMissingET/momentum GenMissingET MissingET
    add Branch UniqueObjectFinder/jets Jet Jet
    add Branch UniqueObjectFinder/electrons Electron Electron
    add Branch UniqueObjectFinder/photons Photon Photon
    add Branch UniqueObjectFinder/muons Muon Muon
    add Branch FatJetFinder/jets FatJet Jet
    add Branch FastJetFinder/jets SmallJet Jet
    add Branch MissingET/momentum MissingET MissingET
    add Branch ScalarHT/energy ScalarHT ScalarHT
    add Branch TrackSmearing/tracks Track Track
    add Branch Rho/rho Rho Rho
    add Branch PileUpMerger/vertices Vertex Vertex
}
