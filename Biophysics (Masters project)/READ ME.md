## Simulating the electrophysiology of the heart

Ectopic beats originating in the pulmonary vein have been found to be a source of atrial fibrillation. Pulmonary vein cells can be categorised into those with pacemaker ability and those
without. Two single cell models of the human pulmonary vein were created for each of these categories by modifying a c++ model of human right atrium. To create the non-pacemaking cell
model, the formulations of the currents for IKr, IKs, IK1, Ito, INa and ICaL were adapted by simulating voltage clamps and modifying the parameters until the output was consistent with experimental
data taken from canine pulmonary vein cells. This model was converted to a pacemaking cell model by adjusting the formulations of ICaL, Ito, INa, IK1 and INCX and adding additional pacemaking
currents If and ICaT that were consistent with data taken from pacemaking rabbit pulmonary vein cells.

### Folder Contents

There is a folder for each of the currents that contribute to the total current propagating through a pulmonary vein cell. These contain c++ and python files that simulate voltage clamp experiments.
Voltage clamp experiments measure current when resting potential is held at a constant value. For each of the currents, the parameters were adjusted until the model fit the experimental data.

These parameters were then added into the 1D tissue model (in the single cell file) to replicate both a pacemaking and non-pacemaking pulmonary vein cell. 
