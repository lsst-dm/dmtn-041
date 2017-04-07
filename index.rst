.. vim: ts=3:sts=3

:tocdepth: 1

.. note::

   This document captures the Winter 2014 design work on the camera geometry (“CamGeom”) system.
   For current documentation, see `the representation of a camera`_ and the `CamGeom documentation`_.

.. _the representation of a camera: https://confluence.lsstcorp.org/display/LSWUG/Representation+of+a+Camera
.. _CamGeom documentation: https://lsst-web.ncsa.illinois.edu/doxygen/x_masterDoxyDoc/afw_camera_geom.html

.. sectnum::

Third party packages
====================

None.

LSST packages used
==================

`afw`_: ``lsst.afw.geom`` (``afwGeom``).

.. _afw: https://github.com/lsst/afw

Requirements
============

#. The implementation will be in Python where possible.
   There will be a (hopefully non-polymorphic) class in C++ to hold the coordinate transforms and other information associated with the detectors and needed in C++.
   This will eliminate any downcasting in Python.
#. The design will define a set of coordinate systems and a way to transform between them.
#. The standard set of coordinates systems will be extensible so that other coordinate systems (and accompanying transforms) can be defined.
#. The design must, at a minimum, support LSST, HSC, Suprime-Cam, MegaCam, SDSS, and DECam.
#. The design will be hierarchical with a top level camera object containing detector objects.
#. The hierarchy will be extensible so that cameras with intermediate levels of sensor grouping (LSST rafts) can iterate over any level.
#. The position and orientation of every sensor may be specified in all 6 axes.
#. Detectors will be iterable from the camera level.
   Every component will be retrievable via index or slot name.
   Component identifiers will contain vendor and serial information.
#. Wavelength-dependent coordinate transforms will be supported, but only minimally, via band-dependent coordianate transforms.
   Monochromatic coordinate transforms (where wavelength is a numeric input) are out of scope because we don't yet understand the use cases well enough to support them.
#. Camera objects will be persistable.
#. The design will provide simple tools for visualizing the camera geometry.
#. The design will make it easy to handle raw data and assembled data and hard to confuse the two.
   We wish to minimize or eliminate the need for flags that indicate if a bbox is in raw or assembled coordinates and whether data has been trimmed.

Coordinate Systems
==================

The following coordinate systems are expected to be supported by every project:

- ``pupil``: x,y radians relative to the optical axis
- ``focalPlane``: x,y mm in the focal plane.
- ``pixels``: nominal x,y unbinned pixels on the entry surface of the detector.

Projects may wish to support additional coordinate systems, such as:

- ``distortedPixels``: x,y unbinned pixel, taking into account such effects as tree rings.
- pupil-like systems at various optical elements (filters, lenses).

.. note::

   Detectors will be positioned with 6 degrees of freedom in the focal plane.
   The actual position of each detector may be reflected in the transformation between pixels and focalPlane.

Classes and API
===============

Note: these are written in pseudocode; details such as pointers (unless essential), references and const are intentionally omitted. All classes are in ``lsst.afw.geom`` or ``lsst.afw.cameraGeom`` unless otherwise noted, and all new classes are in the latter, except ``TransformMap``.

``TransformMap``: C++ (in ``lsst.afw.geom``)
--------------------------------------------

The transform map is a map of ``CoordSys:XYTransform``, where ``CoordSys`` is a template parameter.
Each map will have a native (aka “reference”) coordinate system, and each transform transforms ``Point2D`` objects from the native system to the named system in the forward direction, and supports the reverse transform.
``CoordSys`` will be used as a key in some kind of map (likely ``std::map`` until SWIG can wrap ``unordered_map``), so it will have to follow the appropriate rules for that kind of mapping.

A transform method will transform between any two supported coordinate systems (by transforming the input to the native coordinate system and then transforming that to the output coordinate system).
A variant that can efficiently transform a list of points will also be provided.
To take advantage of caching possibilities, that method will require all input points be in the same coordinate system, so the input is a list of ``Point2D`` instead of ``CameraPoint``.
The output is also a list of a list of ``Point2D``, so a list of points can be transformed multiple times.

If chromatic transforms are desired (e.g. pupil↔︎focalPlane may have a slow color dependence) then the color information must be encoded in the coordinate system.

API
^^^

``TransformMap(CoordSys nativeCoordSys, {CoordSys: XYTransform} transformMap)``
   construct a ``TransformMap``. Adds a unity native→native transform if a native transform is not provided.
``transform(Point2D point, CoordSys fromCoordSys, CoordSys toCoordSys) → Point2D``
   transform a point from one coordinate system to another
``transform(afwGeom.Point2D[] pointList, CoordSys fromCoordSys, CoordSys toCoordSys) → Point2D[]``
   transform a list of points from one coordinate system to another
``getNativeCoordSys() → CoordSys``
   get native coordinate system
``[CoordSys coordSys] → XYTransform``
   get an ``XYTransform`` by coordinate system; the returned transform transforms from native coordinates to the specified ``CoordSys`` in the forward direction
``contains(coordSys) → bool``
   is there an ``XYTransform`` for the specified coordinate system? (Renamed ``__contains__`` in Python so ``coordSys in transformMap`` works)
``getCoordSysList() → CoordSys[]``
   get list of supported coordinate systems
``size() → int``
   get number of transforms (including native->native transform) (renamed ``__len__`` in Python so ``len(transformMap)`` works)

Iteration over the transform mapping will be supported using ``begin``/``end`` in C++ and something suitable in Python.

``AmpInfoCatalog``, ``AmpInfoTable`` and ``AmpInfoRecord``: C++
---------------------------------------------------------------

An afw table of amplifier information, including bounding box, gain, read noise, linearity and sufficient information about raw amplifiers to allow assembling images from raw images. In cases where we use post-ISR data as input (such as our present use of SDSS data) the raw amplifier information may be omitted.

Using an afw table allows users to extend the schema to add additional amplifier-specific information that is needed for a particular camera.

There will be setters and getters for each of the following fields:

Fields
^^^^^^

``name``
   name of amplifier (each amp in a detector must have a unique name, to allow lookup by name)
``bbox``
   bounding box of amp image on assembled image
``gain``
   amplifier gain in e-/ADU
``readNoise``
   amplifier read noise, in e-
``saturation``
   saturation value, in ADU
``readoutCorner``
   readout corner (an enum; one of ``LL``, ``LR``, ``UR``, ``UL``)
``linearityCoeffs``
   linearity coefficients
``linearityType``
   name of linearity algorithm
``hasRawInfo``
   has raw amp information been provided?
``rawBBox``
   bounding box of all pixels on raw image
``rawDataBBox``
   bounding box of data on raw image
``rawFlipX``
   flip X axis when assembling an image?
``rawFlipY``
   flip Y axis when assembling an image?
``rawXYOffset``
   offset for assembling a raw CCD image: desired ``xy0`` - raw ``xy0``.

   This offset is NOT used by ISR; it is primarily for display utilities; it supports construction of a raw CCD image in the case that raw data is provided as individual amplifier images.

   Use ``0,0`` for cameras that supply raw data as a raw CCD image (most cameras).

   Use nonzero for LSST and other cameras that supply raw data as separate amp images with ``xy0=0,0``.
``rawHorizontalOverscanBBox``
   bounding box of valid pixels of horizontal overscan on raw image

   By “valid pixels” we mean to excluded pixels likely to contain electronic artifacts
   and thus make the data unusable for image processing.
``rawVerticalOverscanBBox``
   bounding box of valid pixels of vertical overscan on raw image
``rawPrescanBBox``
   bounding box of valid pixels of (horizontal) prescan on raw image

Here is a pictorial example for ``flipX``, ``flipY``::

    CCD with 4 amps        Desired assembled output     Use these parameters

    --x         x--            y
   |  amp1    amp2 |           |                               flipX       flipY
   y               y           |                       amp1    False       True
                               | CCD image             amp2    True        True
   y               y           |                       amp3    False       False
   |  amp3    amp4 |           |                       amp4    True        False
    --x         x--             ----------- x

This assumes assembled X is always +/- raw X, which is true for CCDs. If some exotic future detector wants to swap X/Y axes then we can add a ``doTranspose`` flag.

``CameraSysPrefix``: C++
------------------------

An incomplete camera coordinate system that only has a coordinate system name (no detector name).
This is used by ``Detector.getCameraSys`` to expand a detector coordinate prefix into a full ``CameraSys`` (at Jim Bosch's excellent suggestion).

A constant will be provided for each detector-specific coordinate system, including ``PIXELS`` and ``ACTUAL_PIXELS`` (non-detector-specific coordinate systems are handled by ``CameraSys``).

API
^^^

``CameraSysPrefix(sysName)``
   construct a ``CameraSysPrefix``
``getSysName() → string``
   coordinate system name, e.g. ``"pixels"``
``operator==(CameraSysPrefix rhs) → bool``
   equality operator
``operator!=(CameraSysPrefix rhs) →bool``
   inequality operator

``CameraSys``: C++
------------------

Class for camera-based coordinate systems; used as the key for ``TransformMap`` in ``Detector`` and ``Camera``.
A constant will be provided for each non-detector-specific coordinate system, including ``FOCAL_PLANE`` and ``PUPIL`` (``CameraSysPrefix`` is used for detector-specific coordinate system prefixes).

API
^^^

``CameraSys(sysName, detectorName="")``
   construct a ``CameraSys``
``getSysName() → string``
   coordinate system name, e.g. "pixels"
``getDetectorName() → string``
   detector name, e.g. "R11 S02"
``hasDetectorName() → bool``
   is detector name non-empty?
``operator==(CameraSys rhs) → bool``
   equality operator
``operator!=(CameraSys rhs) → bool``
   inequality operator
``operator<(CameraSys rhs) → bool``
   to support ``std::map`` (until we switch to ``std::unordered_map``)

``CameraPoint``: C++
--------------------

A struct-ish class that holds:

``point``
   an ``afwGeom.Point2D``
``cameraSys``
   a ``CameraSys``

API
^^^

- ``CameraPoint(afwGeom.Point2D point, CameraSys cameraSys)``
- ``getPoint() → afwGeom.Point2D point``
- ``getCameraSys() → CameraSys cameraSys``

``Orientation``: C++
--------------------

Position, orientation in focal plane.

API
^^^

``Orientation(afwGeom.Point3D fpPosition, afwGeom.Point2D refPoint, afwGeom.Angle yaw, afwGeom.Angle pitch, afwGeom.Angle roll)``
   constructor
``getFpPosition() → afwGeom.Point3D``
   position of detector refPoint in focal plane (mm)
``getRefPoint() → afwGeom.Point2D``
   reference point on detector; offset is measured to this points and all rotations are about this point
``getYaw() → afwGeom.Angle``
   yaw: rotation about :math:`Z` (:math:`X` to :math:`Y`), 1st rotation
``getPitch() → afwGeom.Angle``
   pitch: rotation about :math:`Y'` (:math:`Z'` = :math:`Z` to :math:`X'`), 2nd rotation
``getRoll() → afwGeom.Angle``
   roll: rotation about :math:`X''` (:math:`Y''` = :math:`Y'` to :math:`Z''`), 3rd rotation
``getNQuarter() → int``
   number of quarter rotations about focalPlane Z required to display pixels in a focalPlane mosaic
``makeFpPixelTransform(float pixelSize) → XYTransform``
   make a focalPlane->pixels transform
``makePixelFpTransform(float pixelSize) → XYTransform``
   make a pixels->focalPlane transform

``DetectorType``: C++
---------------------

An enum that identifies the detector type. Possible values are ``SCIENCE``, ``GUIDE``, ``FOCUS`` and ``WAVEFRONT``.

``Detector``: C++
-----------------

``Detector`` holds amplifier information, coordinate transformations important for a detector, and some metadata about the detector including the detector name (which idenitifies the detector slot in the focal plane), detector type, and a serial ID (which idenitifies a particular CCD).

API
^^^

``Detector(string name, int id, DetectorType detectorType, string serial, Box2I bbox, AmpInfoCatalog ampInfo, Orientation orientation, float pixelSize, TransformMap transformMap)``
   constructor
``getId() → int``
   return the detector ID (for use as a key in database tables and such)
``getName() → string``
   return the detector name
``getType() → DetectorType``
   return the detector type
``getSerial() → string``
   serial “number” that identifies the physical detector
``getBBox() → Box2I``
   bounding box of amplifier image (pixels)
``getCorners(CameraSys coordSys) → CameraPoint[4]``
   that describe the extreme corners of the detector in the specified coordinate system
``getCorners(CameraSysPrefix coordSysPrefix) → CameraPoint[4]``
   same for ``CameraSysPrefix``
``transform(CameraPoint cameraPoint, CameraSys toSys) → CameraPoint``
   transform a point to a new coordinate system
``transform(CameraPoint cameraPoint, CameraSysPrefix toSysPrefix) → CameraPoint``
   same for ``CameraSysPrefix``
``getCameraSys(CameraSys cameraSys) → cameraSys``
   return the input unchanged, but check that the detector name matches
``getCameraSys(CameraSysPrefix cameraSysPrefix) → CameraSys``
   return a CameraSys with the detector's name set
``getCenter(CameraSys cameraSys) → CameraPoint``
   get center of detector in specified coordinates
``getCenter(CameraSysPrefix cameraSysPrefix) → CameraPoint``
   same for ``CameraSysPrefix``
``size() → amplifierList.size()`: renamed `__len__``
   in python
``[int index] → AmpInfoRecord``
   get amp info by index
``[string name] → AmpInfoRecord``
   get amp info by amp name
``getAmpInfoCatalog() → AmpInfoCatalog``
   get amp info catalog
``getOrientation() → Orientation``
   get orientation
``getPixelSize() → afwGeom:Point2D``
   x, y pixel size (mm)
``getTransformMap() →  TransformMap``
   return the transform map
``makeCameraPoint(afwGeom:Point2D point, CameraSys cameraSys) → CameraPoint``
   make a ``CameraPoint`` in the given ``CameraSys``
``makeCameraPoint(afwGeom:Point2D point, CameraSysPrefix cameraSysPrefix) → CameraPoint``
   make a ``CameraPoint`` in the given ``CameraSysPrefix``

``DetectorCollection``: Python
------------------------------

A class to hold a collection of ``Detector``\s.
It allows for iteration over all detectors in the collection as well as access by name or index.
Each ``Detector`` must support the same coordinate systems (the constructor will enforce this).

API
^^^

``DetectorCollection(list detectorList)``
   constructor
``__iter__() → iter(detectorList)``
   implement iterator protocol
``__len__() → len(detectorList)``
   implement build in ``len()``
``[int id] → Detector``
   get detector by ID
``[string name] → Detector``
   get detector by name
``getIdIter() → iter(int)``
   iterator over detector IDs
``getNameIter() → iter(str)``
   iterator over detector names
``getCollection(DetectorType detectorType) → DetectorCollection``
   return a collection of detectors of the specified type; this makes it easy to iterate over all science detectors
``getFpBBox() → afwGeom.Box2D``
   get a bounding box that includes all detectors, in focal plane coordinates

``Camera``: Python
------------------

A subclass of ``DetectorCollection`` that also holds a camera ``TransformMap``.
``Camera`` has the ability to transform ``CameraPoints`` to any coordinate system defined in either the camera ``TransformMap`` or its detectors' ``ConversionRegistries``.
``Camera`` will also provide a convenience function to quickly find the detectors that contain a specific ``CameraPoint``.
``Camera`` will also provide class methods to help in building ``Detector``\s.

API
^^^

``Camera(str name, list detectorList, TransformMap transformRegistry)``
   constructor
``getName() → str``
   return the camera name
``findDetectors(CameraPoint cameraPoint) → list of Detectors``
   return a list of ``Detector``\s that contain the given point; typically returns 1 or 0 detectors but may return more if a camera has overlapping detectors
``getTransformMap() → TransformMap``
   return the ``TransformMap`` of the camera
``transform(CameraPoint cameraPoint, CameraSys toCameraSys) → CameraPoint``
   return ``cameraPoint`` transformed to the new coordinate system

   Note that ``toCameraSys`` must have the detector name filled out (if appropriate); ``transform()`` will not search for a detector.
   Call ``findDetectors()`` to search.
``makeCameraPoint(afwGeom:Point2D point, string coordSys) → CameraPoint``
   return a new ``CameraPoint``, raising if ``coordSys`` is not supported by the camera's transform map (if you want a camera point for a particular detector, use ``Detector``)

Possible helper methods for constructing camera/detector objects
================================================================

``ConstructCameraTask``: Python
-------------------------------

This will be a camera specific task that will create a full camera given a time from persisted data.
Following is a Python pseudocode implementation of a ``CameraFactoryTask`` for an arbitrary camera.
There are many ways to read persisted camera data.
This implementation uses a sqlite database containing the detector information.
It implies the following schema:

``Camera``
   - ``Date`` (``DateTime``)
   - ``CameraId`` (``int``)
   - ``pincushion`` (``float``)
   - ``plateScale`` (``float``)

``SlotMap``
   - ``CameraId`` ``int``  (foreign key on ``Camera``)
   - ``DetectorId`` ``int``  (foreign key on ``Detector``)
   - ``SlotName`` ``varchar``
   - ``SlotIndex`` ``int``

``Detector``
   - ``DetectorId`` ``int``
   - ``DetectorSerial`` ``varchar``
   - ``x`` ``float``
   - ``y`` ``float``
   - ``z`` ``float``
   - ``alpha`` ``float``
   - ``beta`` ``float``
   - ``gamma`` ``float``
   - ``pixelSize`` ``float``

``AmpInfo``
   see ``AmpInfoCatalog`` above

Following is code to make a camera given a database like the one above.
Note that the ``makeDetector`` method has been reimplemented in this task.
It could also be defined in the ``Camera`` class as I suggest in the API above.

.. code-block:: py

   import re
   import sqlite3
   import lsst.afw.geom as afwGeom
   from lsst.awf.table import AmpInfoTable, AmpInfoCatalog
   from lsst.afw.cameraGeom import Camera, Detector, DetectorType, Orientation, TransformMap
   import lsst.pex.config as pexConfig
   from datetime import datetime

   class CameraFactoryConfig(pexConfig.Config):
       repoFile = pexConfig.Field(doc="Name of sqlite file to read", dtype=str, default='ccdParams.sqlite')

   class CameraFactoryTask(object):
       configClass = CameraFactoryConfig
       def __init__(self, config):
           self.config = config

           conn = sqlite3.connect(self.config.repoFile)
           self.cur = conn.cursor()

           cameraProps = self.queryCameraProperties(date)
           self.cameraId = cameraProps['cameraId']
           self.pincushion = cameraProps['pincushion']
           self.plateScale = cameraProps['plateScale'] #arcsec/mm

       def run(self, date=datetime.now()):
           '''
           Construct a camera from a database of detector descriptions given a date
           @param date: A datetime object for the time of the camera
           '''
           detectorList = self.queryDetectors()
           pupilTransform = Camera.makePupilFpTransform(self.plateScale, self.pincushion)
           transformRegistry = TransformMap('focalplane', [('pupil', pupilTransform),])
           return Camera(detectorList, transformRegistry)

       def makeDetector(self, ampDictList, detectorType, x, y, z, alpha, beta, gamma,
                        pixelSize, slotName, slotIndex, detectorSerial, *args, **kwargs):
           """
           Make a detector object:
           @param ampDictList: A list of dictionaries, one per amp, that describe the amp properties.
           @detectorType: Science, guiding, ...
           @param x: X Position of detector (mm)
           @param y: Y Position of detector (mm)
           @param z: Z Position of detector (mm)
           @param alpha: First Euler angle of the solid body rotation of the detector
           @param beta: Second Euler angle of the solid body rotation of the detector
           @param gamma: Third Euler angle of the solid body rotation of the detector
           @param pixelSize: The size of pixels for this device (microns)
           @param slotName: The name of the slot
           @param slotIndex: Integer index of the slot
           @param detectorSerial: The serial number of the specific detector.
           """
           ampInfoSchema = AmpInfoTable.makeMinimalSchema()
           ampInfoCatalog = AmpInfoCatalog(ampInfoSchema)

           ccdBox = afwGeom.Box2I()
           for ampDict in ampDictList:
               record = ampInfoCatalog.addNew()
               ccdBox.include(record.getBbox())

           dt = DetectorType(detectorType)

           position = afwGeom.Point3D(x, y, z)
           orientation = Orientation(position, alpha, beta, gamma)
           pixscale = self.plateScale*pixelSize/1000. #convet micron to mm
           fpPixelTransform = Camera.makeFpPixelTransform(orientation, ccdBox, pixscale)
           transformRegistry = TransformMap('pixel', [('focalplane', fpPixelTransform),])
           return Detector(slotName, dt, detectorSerial, ampList, transformRegistry)

       def setAmplifier(self, record, ampBbox, gain, readNoise,...):
           """
           Make an amplifier object.
           @param[in,out] record: record of AmpInfoCatalog to set
           @param ampBbox: A string containing the bounding box of the full amp pixel grid
           @param gain: gain
           ...
           """
           record.setBBox(self.parseBBox(ampBbox, int))
           record.set...

       @staticmethod
       def parseBBox(bboxStr, boxType=int):
           """Parse a bbox string assuming a format like the one used in FITS headers:
              [xLL:yLL,xUR:yUR]
              param boxType: int or float
           """
           typeMap = {int:'I', float:'D'}
           (x0, y0, x1, y1) = re.compile('[\[:,\]]').split(bboxStr[1:-1])
           x0 = boxType(x0)
           y0 = boxType(y0)
           extx = boxType(x1) - x0
           exty = boxType(y1) - y0
           point = getattr(afwGeom, 'Point2%s'%(typeMap[boxType]))
           extent = getattr(afwGeom, 'Extent2%s'%(typeMap[boxType]))
           box = getattr(afwGeom, 'Box2%s'%(typeMap[boxType]))
           return box(point(x0, y0), extent(extx, exty))

       def queryCameraProperties(self, date):
           '''
           Query the database for the camera properties valid for the date.
           @param data: datetime object for the date in question
           @return A dictionary of camera properties
           '''
           dateTempl = "%04i-%02i-%02i %02i:%02i:%02.3f"
           dateStr = dateTempl%(date.year, date.month, date.day, date.hour,
                                date.minute, date.second+date.microsecond/1000000.)
           rows = self.cur.execute('select cameraId, pincushion, plateScale from Camera '+
                                   'where (? - date) > 0 order by (? - date) limit 1',
                                   (dateStr, dateStr))
           row = rows.fetchone()
           return {'cameraId':row[0], 'pincushion':row[1], 'plateScale':row[2]}

       def queryDetectors(self):
           """
           Query a database for the properties of all detectors active at a particular time
           @return detectorList: A list of detector objects
           """

           qstr = 'select s.slotName, s.slotIndex, d.x, d.y, d.z, d.alpha, d.beta, d.gamma, d.serial '+
                  'from Camera c '+
                  'join slotMap s on (c.cameraId = s.cameraId) '+
                  'join Detector d on (s.detectorId = d.detectorId) '+
                  'where c.cameraId = ?'
           rows = self.cur.execute(qstr, (self.cameraId,))

           detectorList = []
           for row in rows:
               detectorProps = dict([(colname[0], row[i]) for i, colname in enumerate(rows.description)])

               amprows = self.cur.execute('select ampName, assembledBbox, gain, readnoise, ampBbox, rawDataBbox, '+
                                  'oscanH, oscanV, flipx, flipy, xOffset, yOffset '+
                                  'from amps where serial = ?', (detectorDict[name]['serial'],))

               ampList = []
               for amprow in amprows:
                   ampList.append(dict([(colname[0], amprow[i]) for i, colname in enumerate(amprows.description)]))
..                detectorList.append(self.makeDetector(ampList, **detectorProps))

           return detectorList

ISR
===

Here is now ISR assembles an image.
Two versions are shown: the first is the default implementation that handles most cameras, and the second handles LSST.

.. code-block:: py

   def assembleImageFromCcdImage(rawImage, ampList, bbox=None):
       """Assemble a CCD image for a camera in which raw data is a single CCD image

       This is how many cameras manage raw data, but not LSSTSim.

       @param[in] rawImage: raw image (Image or MaskedImage). All the amplifier images are included
           somewhere on the raw image. Underscan and overscan make the raw image larger than the assembled image.
       @param[in] ampList: a list of Amplifiers
       @param[in] bbox: bounding box of output amplifier; if None then compute from ampList.
       @return assembled image or maskedImage
       """
       if bbox is None:
           bbox = Box2I()
           for amp in ampList:
               bbox.include(amp.getBBox())

       outImage = Image(bbox)
       for amp in ampList:
           amp.assembleImage(outImage, rawImage)
       return outImage

   def assembleImageFromAmpImages(rawImageList, ampList, bbox=None):
       """Assemble a CCD image for a camera in which raw data is a set of separate amplifier images, e.g. LSSTSim

       @param[in] rawImageList: list of raw images (Image or MaskedImage), one per amplifier
       @param[in] ampList: a list of Amplifiers
       @param[in] bbox: bounding box of output amplifier; if None then compute from ampList.
       @return assembled image or maskedImage
       """
       if len(rawImageList) != len(ampList):
           raise RuntimeError("You must provide one raw image per amplifier")

       if bbox is None:
           bbox = Box2I()
           for amp in ampList:
               bbox.include(amp.getBBox())

       outImage = Image(bbox)
       for rawImage, amp in izip(rawImageList, ampList):
           amp.assembleImage(outImage, rawImage)
       return outImage

Notes
=====

- We assume the coordinate transformations supported by these classes are two dimensional.
  This is not fully general, but we feel it covers the primary use cases and is a very useful simplification.
- ``Camera.transform`` must allow one to specify a particular detector, for example to use in pupil→pixels transforms with a known detector (if not specified then transform calls ``findDetector``.
  ``CameraPoint`` has a detector name field that works well for pixels→pupil, but would require setting to go the other way.
  That may suffice, or we can add a detector argument.
- Our plan is that ``Detector`` can only transform between detector-level coordinate systems, such as ``pixels``, and ``focalPlane`` coordinates; to transform to other camera coordinate systems such as ``pupil`` requires a ``Camera``.
  Supporting additional transforms requires putting a copy of the camera's transform map in each ``Detector`` and slightly more complicated code in ``Detector``'s ``transform`` method to handle the two transform maps.
- Mutability: we plan to make all the above classes immutable after construction.
  This will require one to make a new ``Detector`` or ``Camera`` if either's ``TransformMap`` is to be extended at run time, but we suspect this use case is fairly rare.
- RHL wants enough image assembly and bias subtraction to be in CameraGeom that display code and engineers working on CCDs will not have to use the ISR package.
  Right now there is duplication between ip_isr and CameraGeom and we should take this opportunity to eliminate that.
- RHL wants a minimal “trim and subtract bias” function for image display (cruder than ISR would use).
  That will be written in python.

.. note::

   This document was originally published as an LSST TRAC page at
   https://dev.lsstcorp.org/trac/wiki/Winter2014/Design/CameraGeomDesign.
