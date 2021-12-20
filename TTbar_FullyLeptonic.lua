-- Load the library containing the matrix element
load_modules('/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/TTbar_FullyLeptonic/build/libme_myHappyME.so')

local lepton1 = declare_input("lepton1")
local lepton2 = declare_input("lepton2")
local bjet1 = declare_input("bjet1")
local bjet2 = declare_input("bjet2")

parameters = {
    energy = 13000.,
    top_mass = 173.,
    top_width = 1.491500,
    W_mass = 80.419002,
    W_width = 2.047600,
}

-- Configuration of Cuba
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
    sigma = 0.05,
    -- sigma_range = 5.,
}

lepton1.set_gen_p4("tf_p1::output");

GaussianTransferFunctionOnEnergy.tf_p2 = {
    ps_point = add_dimension(),
    reco_particle = lepton2.reco_p4,
    sigma = 0.05,
    -- sigma_range = 5.,
}

lepton2.set_gen_p4("tf_p2::output")

GaussianTransferFunctionOnEnergy.tf_p3 = {
    ps_point = add_dimension(),
    reco_particle = bjet1.reco_p4,
    sigma = 0.10,
    -- sigma_range = 5.,
}

bjet1.set_gen_p4("tf_p3::output");

GaussianTransferFunctionOnEnergy.tf_p4 = {
    ps_point = add_dimension(),
    reco_particle = bjet2.reco_p4,
    sigma = 0.10,
    -- sigma_range = 5.,
}

bjet2.set_gen_p4("tf_p4::output");

inputs_before_perm = {
    'tf_p1::output',
    'tf_p2::output',
    'tf_p3::output',
    'tf_p4::output',
}

add_gen_permutations(bjet1, bjet2)

inputs = {
  lepton1.gen_p4,
  bjet1.gen_p4,
  lepton2.gen_p4,
  bjet2.gen_p4
}

BreitWignerGenerator.flatter_s13 = {
    ps_point = add_dimension(),
    mass = parameter('W_mass'),
    width = parameter('W_width')
}

BreitWignerGenerator.flatter_s134 = {
    ps_point = add_dimension(),
    mass = parameter('top_mass'),
    width = parameter('top_width')
}

BreitWignerGenerator.flatter_s25 = {
    ps_point = add_dimension(),
    mass = parameter('W_mass'),
    width = parameter('W_width')
}

BreitWignerGenerator.flatter_s256 = {
    ps_point = add_dimension(),
    mass = parameter('top_mass'),
    width = parameter('top_width')
}

BlockD.blockd = {
    p3 = inputs[1],
    p4 = inputs[2],
    p5 = inputs[3],
    p6 = inputs[4],

    pT_is_met = true,

    s13 = 'flatter_s13::s',
    s134 = 'flatter_s134::s',
    s25 = 'flatter_s25::s',
    s256 = 'flatter_s256::s',
}

Looper.looper = {
    solutions = 'blockd::solutions',
    path = Path('boost', 'ttbar', 'integrand')
}

--
-- Start of loop over solutions
--

    full_inputs = {lepton1.gen_p4, 'looper::particles/1', bjet1.gen_p4, lepton2.gen_p4, 'looper::particles/2', bjet2.gen_p4}
    
    BuildInitialState.boost = {
        particles = full_inputs,
        do_transverse_boost = true
    }

    jacobians = { 'flatter_s13::jacobian', 'flatter_s134::jacobian', 'flatter_s25::jacobian', 'flatter_s256::jacobian', 'tf_p1::TF_times_jacobian', 'tf_p2::TF_times_jacobian', 'tf_p3::TF_times_jacobian', 'tf_p4::TF_times_jacobian', 'looper::jacobian' }

    MatrixElement.ttbar = {
      pdf = 'CT10nlo',
      pdf_scale = parameters.top_mass,

      matrix_element = 'myHappyME_sm_P1_Sigma_sm_gg_mupvmbmumvmxbx',
      matrix_element_parameters = {
          card = '/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/TTbar_FullyLeptonic/Cards/param_card.dat'
      },

      initialState = 'boost::partons',

      particles = {
        inputs = full_inputs,
        ids = {
	  {
            pdg_id = -13,
            me_index = 1,
          },
	  {
            pdg_id = 14,
            me_index = 2,
          },
	  {
            pdg_id = 5,
            me_index = 3,
          },
	  {
            pdg_id = 11,
            me_index = 4,
          },
	  {
            pdg_id = -12,
            me_index = 5,
          },
	  {
            pdg_id = -5,
            me_index = 6,
          },
        }
      },

      -- Other jacobians
      jacobians = jacobians
    }

    DoubleLooperSummer.integrand = { input = 'ttbar::output' }

--
-- End of loop over solutions
--

integrand('integrand::sum')

