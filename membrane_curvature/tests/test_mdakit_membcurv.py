"""
Unit and regression test for the membrane_curvature package.
"""

# Import package, test suite, and other packages as needed
import pickle
import pytest
import sys
import mdtraj as md  # This will be gone after refactoring
import itertools as it
from ..lib.mods import core_fast_leaflet, curvature, def_all_beads, mean_curvature, gaussian_curvature, dict2pickle
import numpy as np
from numpy.testing import assert_almost_equal
import MDAnalysis as mda
from membrane_curvature.tests.datafiles import (GRO_MEMBRANE_PROTEIN, XTC_MEMBRANE_PROTEIN, GRO_PO4, XTC_PO4,
                                                GRO_PO4_SMALL, GRO_PO4_MED, GRO_PO4_BIG, GRO_PO4_INVERTED_ID)

# Reference data from datafile
MEMBRANE_CURVATURE_DATA = {
    'mean_curvature': np.array([[0.06396822, 0.2254189, 0.22399141, 0.02498465, 0.05090628,
                                 -0.09999741, -0.10476551, -0.05358351, 0.12631834, 0.09827992],
                                [0.09151143, 0.19973338, 0.13808081, 0.04574428, 0.10289627,
                                 -0.04137634, -0.18304112, -0.2249239, 0.04760685, 0.16729716],
                                [-0.01779114, -0.04597901, -0.00778294, 0.09933345, -0.03402726,
                                 -0.16382893, -0.09622266, -0.15120851, -0.13753756, -0.04029748],
                                [-0.0369991, 0.00233226, 0.16181129, 0.18932471, -0.01349983,
                                 -0.28269057, -0.16712506, 0.01942325, -0.03927295, -0.03811465],
                                [-0.06309625, -0.00820666, 0.16715942, 0.29301601, 0.21669994,
                                 0.07393958, -0.0129063, 0.04319332, 0.14495082, 0.07021294],
                                [-0.05875299, -0.04393517, 0.0837562, -0.05789939, -0.19497179,
                                 0.25517884, 0.25387131, -0.03625653, -0.03502722, -0.01136375],
                                [-0.08988405, -0.10299823, -0.04910499, -0.21748747, -0.41463502,
                                 -0.03913769, 0.1765791, -0.03554145, -0.1677006, -0.10516181],
                                [0.0095539, 0.0364163, 0.00168944, -0.06394463, -0.04456537,
                                 -0.25037463, -0.03814847, -0.02531541, -0.11902046, -0.10177806],
                                [0.00184591, 0.12102329, 0.28902913, 0.0966804, -0.03156109,
                                 -0.16964647, -0.336664, 0.03280685, 0.03212416, -0.08340905],
                                [0.01478105, 0.08136444, 0.23413597, 0.1472945, -0.06672947,
                                 -0.09468121, -0.21140388, 0.03506431, 0.03308529, -0.01943328]]),

    'gaussian_curvature': np.array([[-0.00189129, 0.04934495, 0.04660196, -0.00351411, -0.01943406,
                                     -0.03893078, -0.0246221, -0.03280081, -0.12642188, -0.03525892],
                                    [-0.00693942, 0.03182844, 0.00120618, -0.00352255, 0.00342091,
                                     -0.00142247, 0.02465936, 0.04605395, -0.0106348, -0.00313678],
                                    [-0.01649198, 0.00168855, -0.02731173, -0.01564413, -0.00973393,
                                     0.02368517, -0.00604347, 0.0169452, 0.0113627, -0.00027235],
                                    [0.00088378, -0.00464509, 0.0075489, 0.02877688, -0.00328288,
                                     0.07988983, 0.02487373, -0.00071568, -0.0050775, -0.02188734],
                                    [-0.013283, -0.01107469, 0.02580236, 0.0847512, -0.07205011,
                                     -0.05251712, -0.03977956, -0.03174133, 0.0151017, 0.00278634],
                                    [-0.00688401, -0.01030927, -0.03964658, -0.01066038, 0.01942349,
                                     0.01216024, 0.05635031, -0.0138591, -0.00547026, -0.02372161],
                                    [0.00552771, 0.00111084, -0.0688243, 0.0081551, 0.14126933,
                                     -0.01476609, -0.00715425, -0.0059002, 0.02781192, 0.00485525],
                                    [-0.02954713, -0.01709626, -0.01171343, -0.00766876, -0.01780511,
                                     -0.04009694, -0.00779307, -0.04735893, 0.00799721, -0.00961057],
                                    [-0.0033638, 0.009781, 0.06724205, 0.00795326, 0.00034145,
                                     0.02682387, 0.10108879, -0.01423631, -0.01802192, -0.00922448],
                                    [-0.00313899, -0.00418259, 0.03487913, -0.04456535, -0.00768992,
                                     -0.00642677, 0.0254065, -0.01830984, -0.00904487, -0.01182518]]),

    'z_avg_coords': np.array([[15.52998744, 15.68577027, 15.54460264, 15.42630775, 15.41610875,
                               15.5227999, 15.40030259, 15.47866627, 15.5402739, 15.58982381],
                              [15.7210999, 15.66186774, 15.66852082, 15.5650629, 15.59180416,
                               15.52540419, 15.43453586, 15.31625871, 15.53287837, 15.72669662],
                              [15.65202988, 15.87587967, 15.80443607, 15.2156001, 15.9012679,
                               15.61628637, 15.54388155, 15.51714352, 15.56670715, 15.67794791],
                              [15.79848438, 15.86021965, 15.91864504, 15.69624091, 16.07622673,
                               15.55569959, 15.33016631, 15.51812384, 15.57359306, 15.74745874],
                              [15.99946331, 15.80777305, 15.93547516, 16.06895507, 16.49331617,
                               16.47925099, 15.74151387, 15.63836954, 15.72500532, 15.91917461],
                              [16.00306618, 15.44807954, 14.91708903, 15.62196315, 15.63327226,
                               16.56862831, 16.11996038, 15.89331185, 15.83679604, 15.98079948],
                              [15.76326987, 15.64245421, 15.52017543, 15.47039083, 14.88002689,
                               15.53698847, 15.80894958, 15.80674416, 15.7546208, 15.88995425],
                              [15.70029967, 15.61600719, 15.53771198, 15.47668145, 15.19293903,
                               15.41956353, 15.66099601, 15.71207747, 15.80128901, 15.71542822],
                              [15.68171964, 15.63995413, 15.53009812, 15.21636634, 15.27028742,
                               15.38000841, 15.51385563, 15.64464232, 15.79420718, 15.77794963],
                              [11.9375362, 15.67473274, 15.56693593, 15.28816211, 15.26380259,
                               15.34311786, 15.50101777, 15.5245863, 13.5766693, 15.69553947]]),

    'mean_from_z_avg':
    {'upper': np.array([[0.12555075, 0.02135105, 0.04160992, -0.01648831, 0.00093022,
                         0.05210213, -0.01413144, 0.42855804, 0.11933335, -0.15700576],
                        [0.0119322, -0.03312442, 0.16867713, 0.11186921, -0.06924967,
                         -0.17337403, -0.03075403, 0.00998621, 0.03106606, 0.12265291],
                        [0.01316543, 0.00331553, 0.021467, 0.14910619, -0.02331827,
                         -0.40607687, 0.00254864, 0.09637312, 0.00504627, 0.45395702],
                        [0.07435368, 0.01111333, -0.26369175, 0.02026796, 0.15654274,
                         -0.04698233, -0.07173794, -0.00261281, -0.11980999, -0.1344619],
                        [-0.04629355, -0.00581553, 0.03905765, 0.20953977, 0.3918796,
                         0.0310723, -0.42269884, -0.16012767, -0.01203511, -0.04650698],
                        [-0.00449059, 0.01816231, -0.03467836, -0.09131014, 0.23924521,
                         0.3350186, -0.07902011, -0.16677442, -0.02237877, -0.01612266],
                        [-0.03812502, -0.00993137, -0.05389885, -0.23281703, -0.03940271,
                         0.22679189, 0.15163259, 0.00348903, -0.03725647, 0.23348912],
                        [-0.10594603, -0.14290049, -0.04246863, -0.05403972, -0.14319528,
                         -0.02398386, 0.06424594, 0.03606966, 0.03191477, -0.03725069],
                        [-0.00526536, -0.04075257, -0.04120521, -0.03985792, -0.03383821,
                         -0.03392527, -0.04174206, 0.31247241, 0.19253582, -0.27770588],
                        [0.04900292, 0.03659575, -0.03416543, -0.0529405, 0.00695103,
                         0.01313807, -0.03327879, 0.01243288, 0.04821536, -0.09293237]])},

    'gaussian_from_z_avg':
    {'upper': np.array([[-2.39194800e-02, -3.61634194e-03, -1.08310985e-02,
                         -4.04262352e-02, -1.07953930e-01, -1.58940765e-03,
                         -5.14815931e-02, -4.08411708e-03, -1.73969187e-01,
                         -1.27348350e-02],
                        [-1.80881429e-02, -8.37058009e-04, 1.55266256e-02,
                         7.96503053e-03, -8.70989095e-02, 2.62876632e-02,
                         -5.03604863e-02, -7.13194031e-04, -8.88446427e-01,
                         -1.94381867e-01],
                        [-6.40974462e-03, -9.66907198e-03, 1.30718905e-04,
                         -5.47821599e-02, -1.76480851e-02, 1.54886560e-01,
                         -1.36958801e-02, 2.75305642e-03, -3.89289835e-03,
                         -2.05699114e-02],
                        [-2.13816311e-02, -3.60467810e-03, 6.62439558e-02,
                         -9.62066579e-03, -2.83642243e-02, -6.27955986e-02,
                         -5.31456175e-02, -8.46882808e-03, 1.37664910e-02,
                         1.68963909e-02],
                        [-2.57613213e-03, -5.75152336e-03, -3.04009492e-03,
                         4.38720368e-02, 8.95669062e-02, -2.11138887e-03,
                         1.21192352e-01, 2.35164083e-02, -1.32594642e-02,
                         -8.60185561e-03],
                        [-6.55374523e-03, -7.43313172e-03, -4.29013834e-02,
                         -1.93519324e-03, 2.36387158e-02, 2.59636175e-02,
                         4.16371248e-03, -5.25782227e-03, -2.71611977e-03,
                         3.10811328e-05],
                        [-5.32693655e-03, -2.85760622e-03, -5.49568585e-04,
                         2.04948874e-02, -1.76805356e-02, -1.36782305e-02,
                         -3.16246302e-02, -5.62996471e-03, -1.44329381e-04,
                         -3.73865110e-02],
                        [6.22342127e-03, 1.94289676e-02, 2.82196188e-05,
                         2.69411991e-03, -1.30866347e-02, -2.52293095e-02,
                         -7.22431515e-03, -6.73201764e-03, -2.51286123e-01,
                         -3.22793339e-01],
                        [-2.24613530e-02, 7.82002286e-04, -3.98911117e-04,
                         5.29635968e-04, -4.19699155e-03, -1.30479296e-02,
                         -2.38928651e-04, 4.24992903e-02, 1.70628313e-02,
                         -4.85138439e-02],
                        [-1.94220106e-02, -2.13212988e-04, 1.05435416e-03,
                         1.30710748e-03, -3.49076390e-03, -1.28226749e-02,
                         -1.41168166e-02, -1.08870941e-02, -1.21323346e+00,
                         -1.55778010e-01]])},

    'z_ref': np.array([[15.96339989, 16.36299992, 16.57025003, 16.30183315, 16.14233398,
                        15.94885731, 15.99250078, 16.17940025, 16.08424997, 15.93374944],
                       [15.9995718, 16.48649979, 16.59700012, 16.3276666, 16.26959991,
                        15.77983316, 15.67449999, 15.59950042, 16.05650028, 16.11733341],
                       [15.9301672, 16.04720001, 16.24383338, 16.38975, 16.05666653,
                        15.71950006, 15.7414999, 15.65285724, 15.71783352, 15.91666635],
                       [15.87350019, 16.0994997, 16.45200014, 16.38366667, 15.91100025,
                        15.44099998, 15.55220013, 15.74933386, 15.7957499, 15.93225002],
                       [15.89100003, 16.01559982, 16.45950031, 16.68450022, 16.34674978,
                        16.27950001, 15.97475028, 16.0142498, 16.07933331, 15.96724939],
                       [15.81012511, 15.84583362, 16.09700036, 15.98525, 15.49299908,
                        16.36499977, 16.20639992, 15.8682003, 15.82559967, 15.87400055],
                       [15.73475003, 15.67866707, 15.85220013, 15.60228566, 15.12299967,
                        15.70033328, 15.87920036, 15.80550003, 15.60928576, 15.8010006],
                       [15.79039974, 15.91499996, 15.97549987, 15.80860004, 15.73637486,
                        15.51133362, 15.80240021, 15.78233337, 15.65516663, 15.72000027],
                       [15.8411665, 16.13249969, 16.48759995, 16.25674987, 15.78233369,
                        15.71450011, 15.33062541, 15.99500027, 15.83737516, 15.74500052],
                       [15.82900023, 16.06166649, 16.2973334, 16.43733342, 16.12957178,
                        16.09366608, 15.95349979, 16.22599983, 16.17750025, 16.00225067]]),

    'membrane_protein':
    {'upper':
     {'POPC': [4292, 4304, 4316, 4328, 4340, 4352, 4364, 4376, 4388, 4400, 4412,
               4424, 4436, 4448, 4460, 4472, 4484, 4496, 4508, 4520, 4532, 4544, 4556, 4568,
               4580, 4592, 4604, 4616, 4628, 4640, 4652, 4664, 4676, 4688, 4700, 4712, 4724,
               4736, 4748, 4760, 4772, 4784, 4796, 4808, 4820, 4832, 4844, 4856, 4868, 4880,
               4892, 4904, 4916, 4928, 4940, 4952, 4964, 4976, 4988, 5000, 5012, 5024, 5036,
               5048, 5060, 5072, 5084, 5096, 5108, 5120, 5132, 5144, 5156, 5168, 5180, 5192,
               5204, 5216, 5228, 5240, 5252, 5264, 5276, 5288, 5300, 5312, 5324, 5336, 5348,
               5360, 5372, 5384, 5396, 5408, 5420, 5432, 5444, 5456, 5468, 5480, 5492, 5504,
               5516, 5528, 5540, 5552, 5564, 5576, 5588, 5600, 5612, 5624, 5636, 5648, 5660,
               5672, 5684, 5696, 5708, 5720, 5732, 5744, 5756, 5768, 5780, 5792, 5804, 5816,
               5828, 5840, 5852, 5864, 5876, 5888, 5900, 5912, 5924, 5936, 5948, 5960, 5972,
               5984, 5996, 6008, 6020, 6032, 6044, 6056, 6068, 6080, 6092, 6104, 6116, 6128,
               6140, 6152, 6164, 6176, 6188, 6200, 6212, 6224, 6236, 6248, 6260, 6272, 6284,
               6296, 6308, 6320, 6332, 6344, 6356, 6368, 6380, 6392, 6404, 6416, 6428, 6440,
               6452, 6464, 6476, 6488, 6500, 6512, 6524, 6536, 6548, 6560, 6572, 6584, 6596,
               6608, 6620, 6632, 6644, 6656, 6668, 6680, 6692, 6704, 6716, 6728, 6740, 6752,
               6764, 6776, 6788, 6800, 6812, 6824, 6836, 6848, 6860, 6872, 6884, 6896, 6908,
               6920, 6932, 6944, 6956, 6968, 6980],
      'POPE': [6992, 7004, 7016, 7028, 7040, 7052, 7064, 7076, 7088, 7100, 7112, 7124,
               7136, 7148, 7160, 7172, 7184, 7196, 7208, 7220, 7232, 7244, 7256, 7268, 7280,
               7292, 7304, 7316, 7328, 7340, 7352, 7364, 7376, 7388, 7400, 7412, 7424, 7436,
               7448, 7460, 7472, 7484, 7496, 7508, 7520, 7532, 7544, 7556, 7568, 7580, 7592,
               7604, 7616, 7628, 7640, 7652, 7664, 7676, 7688, 7700, 7712, 7724, 7736, 7748,
               7760, 7772, 7784, 7796, 7808, 7820, 7832, 7844, 7856, 7868, 7880, 7892, 7904,
               7916, 7928, 7940, 7952, 7964, 7976, 7988, 8000, 8012, 8024, 8036, 8048, 8060,
               8072, 8084, 8096, 8108, 8120, 8132, 8144, 8156, 8168, 8180, 8192, 8204, 8216,
               8228, 8240, 8252, 8264, 8276, 8288, 8300, 8312, 8324, 8336, 8348, 8360, 8372,
               8384, 8396, 8408, 8420, 8432, 8444, 8456, 8468, 8480, 8492, 8504, 8516, 8528,
               8540, 8552, 8564, 8576, 8588, 8600, 8612, 8624, 8636, 8648, 8660, 8672, 8684,
               8696, 8708, 8720, 8732, 8744, 8756, 8768, 8780, 8792, 8804, 8816, 8828, 8840,
               8852, 8864, 8876, 8888, 8900, 8912, 8924, 8936, 8948, 8960, 8972, 8984, 8996,
               9008, 9020, 9032, 9044, 9056, 9068, 9080, 9092, 9104, 9116, 9128, 9140, 9152,
               9164, 9176, 9188, 9200, 9212, 9224, 9236, 9248, 9260, 9272, 9284, 9296, 9308,
               9320, 9332, 9344, 9356, 9368, 9380, 9392, 9404, 9416, 9428, 9440, 9452, 9464,
               9476, 9488, 9500, 9512, 9524, 9536, 9548, 9560, 9572, 9584, 9596, 9608, 9620,
               9632, 9644, 9656, 9668, 9680, 9692, 9704]},
     'lower':
     {'POPC': [11532, 11544, 11556, 11568, 11580, 11592, 11604, 11616, 11628, 11640,
               11652, 11664, 11676, 11688, 11700, 11712, 11724, 11736, 11748, 11760,
               11772, 11784, 11796, 11808, 11820, 11832, 11844, 11856, 11868, 11880,
               11892, 11904, 11916, 11928, 11940, 11952, 11964, 11976, 11988, 12000,
               12012, 12024, 12036, 12048, 12060, 12072, 12084, 12096, 12108, 12120,
               12132, 12144, 12156, 12168, 12180, 12192, 12204, 12216, 12228, 12240,
               12252, 12264, 12276, 12288, 12300, 12312, 12324, 12336, 12348, 12360,
               12372, 12384, 12396, 12408, 12420, 12432, 12444, 12456, 12468, 12480,
               12492, 12504, 12516, 12528, 12540, 12552, 12564, 12576, 12588, 12600,
               12612, 12624, 12636, 12648, 12660, 12672, 12684, 12696, 12708, 12720,
               12732, 12744, 12756, 12768, 12780, 12792, 12804, 12816, 12828, 12840,
               12852, 12864, 12876, 12888, 12900, 12912, 12924, 12936, 12948, 12960,
               12972, 12984, 12996, 13008, 13020, 13032, 13044, 13056, 13068, 13080,
               13092, 13104, 13116, 13128, 13140, 13152, 13164, 13176, 13188, 13200,
               13212, 13224, 13236, 13248, 13260, 13272, 13284, 13296, 13308, 13320,
               13332, 13344, 13356, 13368, 13380, 13392, 13404, 13416, 13428, 13440,
               13452, 13464, 13476, 13488, 13500, 13512, 13524, 13536, 13548, 13560,
               13572, 13584, 13596, 13608, 13620, 13632, 13644, 13656, 13668, 13680,
               13692, 13704, 13716, 13728, 13740, 13752, 13764, 13776, 13788, 13800,
               13812, 13824, 13836, 13848, 13860, 13872, 13884, 13896, 13908, 13920,
               13932, 13944, 13956, 13968, 13980, 13992, 14004, 14016, 14028, 14040,
               14052, 14064, 14076, 14088, 14100, 14112, 14124, 14136, 14148, 14160,
               14172, 14184, 14196, 14208, 14220, 14232, 14244, 14256, 14268, 14280],

      'POPE': [14292, 14304, 14316, 14328, 14340, 14352, 14364, 14376, 14388, 14400, 14412,
               14424, 14436, 14448, 14460, 14472, 14484, 14496, 14508, 14520, 14532, 14544,
               14556, 14568, 14580, 14592, 14604, 14616, 14628, 14640, 14652, 14664, 14676,
               14688, 14700, 14712, 14724, 14736, 14748, 14760, 14772, 14784, 14796, 14808,
               14820, 14832, 14844, 14856, 14868, 14880, 14892, 14904, 14916, 14928, 14940,
               14952, 14964, 14976, 14988, 15000, 15012, 15024, 15036, 15048, 15060, 15072,
               15084, 15096, 15108, 15120, 15132, 15144, 15156, 15168, 15180, 15192, 15204,
               15216, 15228, 15240, 15252, 15264, 15276, 15288, 15300, 15312, 15324, 15336,
               15348, 15360, 15372, 15384, 15396, 15408, 15420, 15432, 15444, 15456, 15468,
               15480, 15492, 15504, 15516, 15528, 15540, 15552, 15564, 15576, 15588, 15600,
               15612, 15624, 15636, 15648, 15660, 15672, 15684, 15696, 15708, 15720, 15732,
               15744, 15756, 15768, 15780, 15792, 15804, 15816, 15828, 15840, 15852, 15864,
               15876, 15888, 15900, 15912, 15924, 15936, 15948, 15960, 15972, 15984, 15996,
               16008, 16020, 16032, 16044, 16056, 16068, 16080, 16092, 16104, 16116, 16128,
               16140, 16152, 16164, 16176, 16188, 16200, 16212, 16224, 16236, 16248, 16260,
               16272, 16284, 16296, 16308, 16320, 16332, 16344, 16356, 16368, 16380, 16392,
               16404, 16416, 16428, 16440, 16452, 16464, 16476, 16488, 16500, 16512, 16524,
               16536, 16548, 16560, 16572, 16584, 16596, 16608, 16620, 16632, 16644, 16656,
               16668, 16680, 16692, 16704, 16716, 16728, 16740, 16752, 16764, 16776, 16788,
               16800, 16812, 16824, 16836, 16848, 16860, 16872, 16884, 16896, 16908, 16920,
               16932, 16944, 16956, 16968, 16980, 16992, 17004, 17016, 17028, 17040]}},

    'small': {'lower': {'POPC': [5, 6, 7], 'POPE': [8, 9]}, 'upper': {'POPC': [0, 1, 4], 'POPE': [2, 3]}},

    'med': {'lower': {'POPC': [12, 13, 14, 15, 16],
                      'POPE': [17, 18, 19, 20, 21, 22, 23, 24]},
            'upper': {'POPC': [0, 1, 2, 3, 4, 5, 6, 7, 11],
                      'POPE': [8, 9, 10]}},

    'big': {'lower': {'POPC': [25, 26, 27, 28, 29, 30, 31, 32, 33, 34],
                      'POPE': [35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49]},
            'upper': {'POPC': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 23, 24],
                      'POPE': [15, 16, 17, 18, 19, 20, 21, 22]}}
}


@ pytest.fixture()
def universe():
    u = mda.Universe(GRO_MEMBRANE_PROTEIN, XTC_MEMBRANE_PROTEIN)
    return u


@ pytest.fixture()
def mdtraj_po4():
    # trajectory PO4 beads only imported using mdtraj, gone after refactoring
    mdtraj = md.load(XTC_PO4, top=GRO_PO4)
    return mdtraj


@ pytest.fixture()
def md_ref_beads():
    # reference beads using mdtraj, gone after refactoring. Select upper leaflet
    topology = md.load(GRO_PO4).topology
    md_ref_beads = {'upper': {'POPC': topology.select('resname POPC and index 1 to 457').astype(int).tolist()}}
    return md_ref_beads


def test_membrane_curvature_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "membrane_curvature" in sys.modules


def test_dict_to_pickle(tmpdir):
    name = 'test_pickle_output'
    dict_ = {'A': 1, 'B': 2, 'C': 3}
    with tmpdir.as_cwd():
        dict2pickle(name, dict_)
        unpickled = pickle.load(open(name + '.pickle', 'rb'))
        assert dict_ == unpickled


def test_gaussian_curvature():
    K_test = gaussian_curvature(MEMBRANE_CURVATURE_DATA['z_ref'])
    for k, k_test in zip(MEMBRANE_CURVATURE_DATA['gaussian_curvature'], K_test):
        assert_almost_equal(k, k_test)


def test_mean_curvature():
    H_test = mean_curvature(MEMBRANE_CURVATURE_DATA['z_ref'])
    for h, h_test in zip(MEMBRANE_CURVATURE_DATA['mean_curvature'], H_test):
        assert_almost_equal(h, h_test)


def test_core_fast_leaflet(md_ref_beads, mdtraj_po4):
    jump = 1
    n_cells = 10
    max_width = 19
    z_calc = np.zeros([n_cells, n_cells])
    core_fast_leaflet(z_calc, "upper", mdtraj_po4, jump, n_cells, ["POPC"], md_ref_beads, max_width)
    for z, z_test in zip(MEMBRANE_CURVATURE_DATA['z_avg_coords'], z_calc):
        assert_almost_equal(z, z_test)


def test_curvature(md_ref_beads, mdtraj_po4):
    jump = 1
    n_cells = 10
    max_width = 19
    z_test = np.zeros([n_cells, n_cells])
    lf = 'upper'
    core_fast_leaflet(z_test, lf, mdtraj_po4, jump, n_cells, ["POPC"], md_ref_beads, max_width)
    z_ref = {lf: z_test}
    H_test, K_test = curvature(z_ref, [lf], n_cells)

    for h, h_test in zip(MEMBRANE_CURVATURE_DATA['mean_from_z_avg'][lf], H_test[lf]):
        assert_almost_equal(h, h_test)

    for k, k_test in zip(MEMBRANE_CURVATURE_DATA['gaussian_from_z_avg'][lf], K_test[lf]):
        assert_almost_equal(k, k_test)


def test_def_all_beads_membrane_protein():
    lipid_types = ['POPC', 'POPE']
    lfs = ["upper", "lower"]
    head_list = [4286, 11531, 18891]
    topology = md.load(GRO_MEMBRANE_PROTEIN).topology
    beads_test = def_all_beads(lipid_types, lfs, head_list, topology)
    for lt, lf in zip(lipid_types, lfs):
        for bead, bead_t in zip(MEMBRANE_CURVATURE_DATA['membrane_protein'][lf][lt], beads_test[lf][lt]):
            assert int(bead) == int(bead_t)


def test_def_all_beads_small_gro():
    lipid_types = ['POPC', 'POPE']
    lfs = ["upper", "lower"]
    head_list = [0, 4, 11]
    topology = md.load(GRO_PO4_SMALL).topology
    beads_test = def_all_beads(lipid_types, lfs, head_list, topology)
    for lt, lf in zip(lipid_types, lfs):
        for bead, bead_t in zip(MEMBRANE_CURVATURE_DATA['small'][lf][lt], beads_test[lf][lt]):
            assert int(bead) == int(bead_t)


def test_def_all_beads_med_gro():
    lipid_types = ['POPC', 'POPE']
    lfs = ["upper", "lower"]
    head_list = [0, 11, 25]
    topology = md.load(GRO_PO4_MED).topology
    beads_test = def_all_beads(lipid_types, lfs, head_list, topology)
    for lt, lf in zip(lipid_types, lfs):
        for bead, bead_t in zip(MEMBRANE_CURVATURE_DATA['med'][lf][lt], beads_test[lf][lt]):
            assert int(bead) == int(bead_t)


def test_def_all_beads_big_gro():
    lipid_types = ['POPC', 'POPE']
    lfs = ["upper", "lower"]
    head_list = [0, 24, 50]
    topology = md.load(GRO_PO4_BIG).topology
    beads_test = def_all_beads(lipid_types, lfs, head_list, topology)
    for lt, lf in zip(lipid_types, lfs):
        for bead, bead_t in zip(MEMBRANE_CURVATURE_DATA['big'][lf][lt], beads_test[lf][lt]):
            assert int(bead) == int(bead_t)


@pytest.mark.xfail(reason='Non sequential indexes in leaflets. Gone after enhancement of atom selection.')
def test_non_linear_indexes():
    lipid_types = ['POPC', 'POPE']
    lfs = ["upper", "lower"]
    head_list = [0, 4, 11]
    topology = md.load(GRO_PO4_INVERTED_ID).topology
    beads_test = def_all_beads(lipid_types, lfs, head_list, topology)
    for lt, lf in zip(lipid_types, lfs):
        for bead, bead_t in zip(MEMBRANE_CURVATURE_DATA['small'][lf][lt], beads_test[lf][lt]):
            assert int(bead) == int(bead_t)
