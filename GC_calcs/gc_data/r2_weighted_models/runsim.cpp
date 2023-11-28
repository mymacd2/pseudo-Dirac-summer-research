
#include "CRPropa.h"

#include "ppInteractions/ParticleDecay.h"
#include "ppInteractions/NucleusNucleusInteraction.h"

#include "observer_utils.h"
#include "gasmaps.h"

#include <iostream>

using namespace crpropa;

ref_ptr<Candidate> Source::getCandidate() const {

        double R_min = 0.3 * kpc; 

	ref_ptr<Candidate> candidate = new Candidate();
        double w = 1.0; 
        bool reject_candidate = true;
        while (reject_candidate) {

                for (int i = 0; i < features.size(); i++)
                        (*features[i]).prepareCandidate(*candidate);

                // reject_candidate = false;
                Vector3d sun_pos = Vector3d(8.5*kpc, 0., 0.);
                Vector3d pos = candidate->source.getPosition();
                double R = (sun_pos - pos).getR();

                if (R < R_min)   
                        reject_candidate = false;
                else if ( Random::instance().rand() * pow(R_min, -2) <= pow(R, -2) ) {
                        reject_candidate = false;
                        w = pow( R/R_min, 2);
                }                        
        }
        candidate->setWeight(w);
	return candidate;
}

int main(int argc, char* argv[]) {

        // command line arguments: number of particles to inject, output_file 
        if (argc < 3) { return 1; }
        std::string output_file = argv[1];
        int N = static_cast<int>(std::stod(argv[2]));

        std::cout << "\nsaving outputs to: " << output_file << std::endl;
        std::cout << "injecting " << N << " particles!\n" << std::endl;

        double thinning = 0.;           
        ModuleList sim;

        // magnetic field
        // double Bphi = 3. * microgauss; 
        // ref_ptr<Cylindrical_UniformMagneticField> B = new Cylindrical_UniformMagneticField(0., Bphi, 0.);
        ref_ptr<JF12Field> B = new JF12Field(); 

        // gas distribution
        // double gas = 1.0 * pow(1/cm, 3); 
        ref_ptr<Grid1f> gas = get_GALPROP_gasmap();

        // propagation
        double tol = 1e-4;
        double minstep = 1*pc;
        double maxstep = 1*kpc;
        // ref_ptr<PropagationBP> prop = new PropagationBP( B, tol, minstep, maxstep );

        double epsilon = 0.1; 
        ref_ptr<DiffusionSDE> prop = new DiffusionSDE( B, tol, minstep, maxstep, epsilon );
        // prop->setAlpha(1./3.);
        // prop->setScale(1.);
        sim.add( prop );

        // interactions
        ref_ptr<NucleusNucleusInteraction> ppint = new NucleusNucleusInteraction(gas);
        ppint->setThinning(thinning);
        sim.add( ppint );

        ref_ptr<ParticleDecay> pdint = new ParticleDecay();
        pdint->setThinning(thinning);
        sim.add( pdint );

        // conditions
        ref_ptr<MinimumEnergyPerParticleId> minimum_E = new MinimumEnergyPerParticleId();
        minimum_E->add( nucleusId(1,1), 10*GeV );
        // deactivate electrons and photons immediately to save time. 
        minimum_E->add(  11, 1e3*EeV );
        minimum_E->add( -11, 1e3*EeV );
        minimum_E->add(  22, 1e3*EeV );        
        sim.add( minimum_E );

        double rMax = 20 * kpc;
        double zMax = 2 * kpc; 
        Vector3d origin = Vector3d(0.);
        sim.add( new CylindricalBoundary( origin, 2*zMax, rMax) );

        // output: 
        ref_ptr<TextOutput> output = new TextOutput(output_file);
        output->setLengthScale(kpc);
        output->setEnergyScale(eV);

        output->disableAll();
        output->set(Output::OutputColumn::TrajectoryLengthColumn, true);
        output->set(Output::OutputColumn::CurrentIdColumn, true);
        output->set(Output::OutputColumn::CurrentEnergyColumn, true);
        output->set(Output::OutputColumn::CurrentPositionColumn, true);
        output->set(Output::OutputColumn::SourceEnergyColumn, true);                    // for reweighting to PL
        // output->set(Output::OutputColumn::SourcePositionColumn, true);
        // output->set(Output::OutputColumn::SerialNumberColumn, true);
        output->set(Output::OutputColumn::WeightColumn, true);
        
        // observer for neutrinos: 
        ref_ptr<Observer> obs = new Observer();
        // obs->add( new ObserverDetectNucleus() ); 
        obs->add( new ObserverDetectNeutrino() ); 
        obs->add( new ObserverInactiveVeto() );
        obs->setDeactivateOnDetection( true );
        obs->onDetection(output);
        sim.add(obs);

        // source 
        ref_ptr<Source> source = new Source();
        source->add( new SourceParticleType( nucleusId(1,1) ) );
        source->add( new SourcePowerLawSpectrum(TeV, EeV, -1)) ;
        source->add( new SourceIsotropicEmission() );
        source->add( new SourcePulsarDistribution() );

        sim.setShowProgress(true);

        sim.run(source, N, true, false);    // recursive, secondariesFirst
        output->close();
        
        return 0;
}