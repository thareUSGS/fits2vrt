# fits2vrt
WCS FITS conversion to Virtual GDAL header

### Howto
From your python prompt:

```python
>>> import fits2vrt
>>> a=fits2vrt.fitskeys('myimage.fits')
>>> a.fits2vrt()
```

will create `myimage.vrt`

### Dependencies
To have `fits2vrt` working you need to install
* [Astropy](http://www.astropy.org/)
* [GDAL for python](https://pypi.python.org/pypi/GDAL/)
