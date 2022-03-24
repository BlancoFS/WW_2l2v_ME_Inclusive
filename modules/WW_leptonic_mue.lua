-- Load the library containing the matrix element
load_modules('/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/WW_leptonic_noJet_ME/build/libme_WW_leptonic_noJet_ME.so')

local lepton1 = declare_input("lepton1")
local lepton2 = declare_input("lepton2")

parameters = {
    energy = 13000.,
    top_mass = 173.,
    top_width = 1.491500,
    W_mass = 80.419002,
    W_width = 2.047600
}

cuba = {
    relative_accuracy = 0.05,
    verbosity = 0,
    max_eval = 10000,
    n_start = 500,
    -- algorithm = "divonne"
}

GaussianTransferFunctionOnEnergy.tf_p1 = {
    ps_point = add_dimension(),
    reco_particle = lepton1.reco_p4,
    sigma = 0.10,
    -- sigma_range = 5.,
}

lepton1.set_gen_p4("tf_p1::output");

GaussianTransferFunctionOnEnergy.tf_p2 = {
    ps_point = add_dimension(),
    reco_particle = lepton2.reco_p4,
    sigma = 0.10,
    -- sigma_range = 5.,
}

lepton2.set_gen_p4("tf_p2::output");

BreitWignerGenerator.flatter_s13 = {
    ps_point = add_dimension(),
    mass = parameter('W_mass'),
    width = parameter('W_width')
}

BreitWignerGenerator.flatter_s24 = {
    ps_point = add_dimension(),
    mass = parameter('W_mass'),
    width = parameter('W_width')
}

BlockF.blockf = {
    p3 = lepton1.reco_p4,
    p4 = lepton2.reco_p4;

    s13 = 'flatter_s13::s',
    s24 = 'flatter_s24::s',
    q1 = add_dimension(),
    q2 = add_dimension()
}

Looper.looper = {
    solutions = 'blockf::solutions',
    path = Path('initial_state', 'WW', 'integrand')
}

--
-- Start of loop over solutions
--

    full_inputs = { lepton2.reco_p4, 'looper::particles/2', lepton1.reco_p4, 'looper::particles/1' }
    
    BuildInitialState.initial_state = {
        particles = full_inputs,
        do_transverse_boost = true
    }

    jacobians = { 'flatter_s13::jacobian', 'flatter_s24::jacobian', 'tf_p1::TF_times_jacobian', 'tf_p2::TF_times_jacobian', 'looper::jacobian' }

    MatrixElement.WW = {
      pdf = 'CT10nlo',
      pdf_scale = parameter('W_mass'),

      matrix_element = 'WW_leptonic_noJet_ME_sm_P1_Sigma_sm_uux_epvemumvmx',
      matrix_element_parameters = {
          card = '/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/WW_leptonic_noJet_ME/Cards/param_card.dat'
      },

      initialState = 'initial_state::partons',

      particles = {
        inputs = full_inputs,
        ids = {
          {
            pdg_id = -11,
            me_index = 1,
          },
          {
            pdg_id = 12,
            me_index = 2,
          },
          {
            pdg_id = 13,
            me_index = 3,
          },
          {
            pdg_id = -14,
            me_index = 4,
          },
        }
      },

      -- Other jacobians
      jacobians = jacobians
    }

    DoubleLooperSummer.integrand = { input = 'WW::output' }

integrand('integrand::sum')