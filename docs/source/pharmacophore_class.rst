Pharmacophore class
===================

.. currentmodule:: pmapper.pharmacophore


.. rubric:: Class

.. autosummary::
   :nosignatures:
   
   Pharmacophore


.. rubric:: Initialization and setup

.. autosummary::
   :nosignatures:

   Pharmacophore.__init__
   Pharmacophore.update   


.. rubric:: Create pharmacophore

.. autosummary::
   :nosignatures:
   
   Pharmacophore.load_from_mol
   Pharmacophore.load_from_smarts
   Pharmacophore.load_from_feature_factory
   Pharmacophore.load_from_atom_ids


.. rubric:: Pharmacophore properties and methods

.. autosummary::
   :nosignatures:

   Pharmacophore.get_bin_step
   Pharmacophore.get_descriptors
   Pharmacophore.get_feature_coords
   Pharmacophore.get_features_count
   Pharmacophore.get_fp
   Pharmacophore.get_fp2
   Pharmacophore.get_graph
   Pharmacophore.get_mirror_pharmacophore
   Pharmacophore.get_mol
   Pharmacophore.get_signature_md5
   Pharmacophore.iterate_pharm
   Pharmacophore.iterate_pharm1


.. rubric:: Pharmacophore match

.. autosummary::
   :nosignatures:

   Pharmacophore.fit_model
   

.. rubric:: Load/save file

.. autosummary::
   :nosignatures:

   Pharmacophore.load_from_pma
   Pharmacophore.load_from_xyz
   Pharmacophore.load_ls_model
   Pharmacophore.load_from_file
   Pharmacophore.save_ls_model
   Pharmacophore.save_to_pma
   Pharmacophore.save_to_xyz
   

.. autoclass:: Pharmacophore
   :members:
   :special-members: __init__
   :inherited-members:
