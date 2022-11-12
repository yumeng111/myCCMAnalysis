.. _recclasses:

Recclasses
==========
This project is a self-contained collection of storage classes inserted into frames by modules other projects, including their pybindings and table converters.

**Maintainer:** Kevin Meagher

.. toctree::
   :maxdepth: 1

   parameter_storage
   laputop_params
   release_notes

Overview
--------
This project hosts a collection of storage classes, derived from I3FrameObject. These classes are filled by other projects (for example, I3LaputopParams is produced by project toprec). It was decided to put all these storage classes into this self-contained project with minimal dependencies, to offer convenient read/write access for consumers, who just want to do data analysis based on this additional information or who just want to browse events with Steamshovel. In that sense, recclasses is the analog of simclasses.

A noteworthy class in this project is the :ref:`parameter_storage`, a general storage solution for users who want to store fitted parameters with full covariance matrix. Currently, it is used by :ref:`laputop_params`.

Most storage classes also have converter. To learn how to convert content from I3 files into table format, read :ref:`hdfwriter` or :ref:`rootwriter`.

.. _`hdfwriter`: http://software.icecube.wisc.edu/offline_trunk/projects/hdfwriter/
.. _`rootwriter`: http://software.icecube.wisc.edu/offline_trunk/projects/rootwriter/

Storage classes
---------------

Also see the full :ref:`doxygen <recclasses>` docs.

===========================================  ==================================================  =========================
storage classes                              converter                                           associated project
===========================================  ==================================================  =========================
:cpp:class:`I3CLastFitParams`                            I3CLastFitParamsConverter               :ref:`clast`
:cpp:class:`CramerRaoParams`                              CramerRaoParamsConverter               :ref:`cramer-rao`
:cpp:class:`I3CscdLlhFitParams`                        I3CscdLlhFitParamsConverter               :ref:`cscd-llh`
:cpp:class:`I3DST`                                                                               :ref:`dst`
:cpp:class:`I3DST13`                                                                             :ref:`dst`
:cpp:class:`I3DST16`                                                                             :ref:`dst`
:cpp:class:`I3DSTHeader`                                                                         :ref:`dst`
:cpp:class:`I3DSTHeader13`                                                                       :ref:`dst`
:cpp:class:`I3DSTHeader16`                                                                       :ref:`dst`
:cpp:class:`I3DSTReco13`                                                                         :ref:`dst`
:cpp:class:`I3DipoleFitParams`                          I3DipoleFitParamsConverter               :ref:`dipolefit`
:cpp:class:`I3DirectHitsValues`                        I3DirectHitsValuesConverter (Python)      :ref:`CommonVariables`
:cpp:class:`I3FillRatioInfo`                              I3FillRatioInfoConverter               :ref:`fill-ratio`
:cpp:class:`I3FiniteCuts`                                                                        :ref:`finiteReco`
:cpp:class:`I3HitMultiplicityValues`              I3HitMultiplicityValuesConverter (Python)      :ref:`CommonVariables`
:cpp:class:`I3HitStatisticsValues`                  I3HitStatisticsValuesConverter (Python)      :ref:`CommonVariables`
:cpp:class:`I3LaputopParams`                              I3LaputopParamsConverter               :ref:`toprec`
:cpp:class:`I3LineFitParams`                              I3LineFitParamsConverter               :ref:`linefit`
:cpp:class:`OMKeyLink`
:cpp:class:`OMKeyLinkSet`
:cpp:class:`OMKeyPair`
:cpp:class:`OMKeySet`
:cpp:class:`I3OpheliaFirstGuessTrack`            I3OpheliaFirstGuessTrackConverter               :ref:`Ophelia`
:cpp:class:`I3OpheliaParticle`                                                                   :ref:`Ophelia`
:cpp:class:`I3OpheliaRecoResult`                                                                 :ref:`Ophelia`
:cpp:class:`I3ParticleIntersections`                                                             :ref:`VHESelfVeto`
:cpp:class:`I3PortiaEvent`                                  I3PortiaEventConverter               :ref:`Portia`
:cpp:class:`I3PortiaPulse`                                                                       :ref:`Portia`
:cpp:class:`I3STConfiguration`                                                                   :ref:`STTools`
:cpp:class:`I3ShieldHitRecord`                          I3ShieldHitRecordConverter               :ref:`shield`
:cpp:class:`I3StartStopParams`                                                                   :ref:`finiteReco`
:cpp:class:`I3TensorOfInertiaFitParams`        I3TensorOfInertiaFitParamsConverter               :ref:`tensor-of-inertia`
:cpp:class:`I3TimeCharacteristicsValues`      I3TimeCharacteristicsValuesConverter (Python)      :ref:`CommonVariables`
:cpp:class:`I3TopLateralFitParams`                  I3TopLateralFitParamsConverter               :ref:`toprec`
:cpp:class:`I3TopRecoPlaneFitParams`              I3TopRecoPlaneFitParamsConverter               :ref:`toprec`
:cpp:class:`I3TrackCharacteristicsValues`    I3TrackCharacteristicsValuesConverter (Python)      :ref:`CommonVariables`
:cpp:class:`I3Veto`                                                I3VetoConverter (Python)      :ref:`CascadeVariables`
:cpp:class:`I3VetoShort`                                                                         :ref:`CascadeVariables`
===========================================  ==================================================  =========================
