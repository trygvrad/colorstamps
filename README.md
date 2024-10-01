# colorstamps

Python package for 2D colormaps. 
Most included colormaps are based on the 'CAM02-LCD' colorspace as defined in the package colorspacious (https://pypi.org/project/colorspacious/)

Documentation is hosted at https://colorstamps.readthedocs.io/en/latest/ but a summary follows below:


Installation:
pip install colorstamps

Example:
```python
import matplotlib.pyplot as plt
import colorstamps

img = colorstamps.helpers.get_random_data() # numpy array of shape (100,200,2) with 2d data to plot    
rgb, stamp = colorstamps.apply_stamp(img[:,:,0], img[:,:,1], 'peak',
                                   vmin_0 = -1.2, vmax_0 = 1.2,
                                   vmin_1 = -1, vmax_1 = 1,
                                 )

fig, axes = plt.subplots(1,2,figsize=(10,3), dpi = 100)    
axes[0].imshow(rgb)

# show colormap as overlay
overlaid_ax = stamp.overlay_ax(axes[0], lower_left_corner = [0.7,0.85], width = 0.2)
overlaid_ax.set_ylabel(r'$\phi$')
overlaid_ax.set_xlabel(r'$\omega$')

# also show colormap as in separate ax to illustrate functionality
stamp.show_in_ax(axes[1])
axes[1].set_ylabel(r'$\phi$')
axes[1].set_xlabel(r'$\omega$')
```


![](docs/source/images/example0.png?raw=true)

Examples of included colormaps are shown below:

![](docs/source/images/colormaps.png?raw=true)
### Customization
All the colormaps are called by colorstamps.stamps.get_cmap() and the following keyword can be used with either colorstamps.apply_stamp() or colorstamps.stamps.get_cmap() to customize the colormaps:
```
l: int, the size of the colormap will be (l,l), defaults to 256 if None
  rot: float, rotation of the colormap (where applicable)
  J: array-like of length 2 (float,float), determins min and max luminocity where applicable
  sat: float, maximum saturation where applicable
  limit_sat: string, 'individual' or 'shared'. How saturation is limited for relevant colormaps when colors outside sRGB are required
              'individual': each combination J, hue in the colormap has an individual limit to saturation
              'shared': for each J, all hues share a limit, the maximum where all hues can be represented
  a: range along a-axis, array-like [min,max]
              Used to move the center of the colormap where applicable.
              Defaults to (-1,1) which is then multiplied internally with sat
  b: range along b-axis, see a.
Additionally, for radial colormaps the name may have a postfix separated by a space, i.e. 'cone tr'
  the postfix must be some combination of (t,b) and/or (l,r) which defines the quadrant/side of the colormap to include
  t-> top, b-> bottom, r-> right, l-> left, and 'tr'-> top right, etc.
```

### Evaluation of colormaps

The package also includes a method for evaluating colormaps, by evaluating the distribution of lighness, hue, saturation, 
and how it can be percieved by those with partial colorblindness
```python
stamp = colorstamps.Stamp('hsv')
fig, ax = stamp.eval()
stamp = colorstamps.Stamp('peak')
fig, ax = stamp.eval()
```
![](docs/source/images/eval_hsv.png?raw=true)
![](docs/source/images/eval_peak.png?raw=true)

The package also supports different methods for clipping data to the colormap, using the 'clip' keyword in colorstamps.apply_stamp()

![](docs/source/images/point_outside_colormap.png?raw=true)

### Additional colormaps

Additional colormaps available by converting 1d colormaps in matplotlib to 2d colormaps by varying the lightness along the y-axis:

![](docs/source/images/mpl_colormaps.png?raw=true)

### User-provided colormaps
Custom colormaps may be integrated by providing a numpy array of shape (l,l,3) detailing a 2d colormap instead of a name when calling  colorstamps.apply_stamp()
```python
my_cmap = np.zeros((256,256,3))
my_cmap[:,:,0] = np.linspace(0,1,256)[:,np.newaxis]
my_cmap[:,:,2] = np.linspace(0,1,256)[np.newaxis,:]
my_cmap[:,:,1] = 0.5*(my_cmap[:,:,0]+my_cmap[:,:,2])

img = colorstamps.helpers.get_random_data() # numpy array of shape (100,200,2) with 2d data to plot    
fig, ax = plt.subplots(1,1,figsize=(5,3), dpi = 100)    
rgb, stamp = colorstamps.apply_stamp(img[:,:,0], img[:,:,1], my_cmap)

ax.imshow(rgb)
overlaid_ax = stamp.overlay_ax(ax, lower_left_corner = [0.66,0.85], width = 0.2)
overlaid_ax.set_ylabel(r'$\phi$')
overlaid_ax.set_xlabel(r'$\omega$')
```
![](docs/source/images/custom_cmap.png?raw=true)

### Colormaps for line plots
Colormaps can also be used in combination with line plots:
```python
# intensities: array of shape (l,l,n)
# omega: array of shape (l,l)
# phi: array of shape (l,l)
# phase: array of shape (n)

# here l = 5, and n is 41

var_0 = np.array(omega)
var_1 = np.array(phi)
l = phi.shape[0]
rgb, stamp = colorstamps.apply_stamp(var_0, var_1, 'abyss', l = l,
	   vmin_0 = np.min(var_0) - 0.5*(np.max(var_0)-np.min(var_0))/l, 
	   vmax_0 = np.max(var_0) + 0.5*(np.max(var_0)-np.min(var_0))/l,
	   vmin_1 = np.min(var_1) - 0.5*(np.max(var_1)-np.min(var_1))/l, 
	   vmax_1 = np.max(var_1) + 0.5*(np.max(var_1)-np.min(var_1))/l,
	 )


fig, ax = plt.subplots(dpi = 300)
for i, row in enumerate(intensities):
    for j, intensity in enumerate(row):
	    ax.plot(phase, intensity/np.max(intensity), color = rgb[i,j])
	    
ax.set_xticks(phase[::10])
ax.set_xticklabels(phase[::10]%360)
ax.set_ylabel('Amplitude')
ax.set_xlabel('Phase')
ax.set_xlim((0,360))
overlaid_ax = stamp.overlay_ax(ax, lower_left_corner = [0.15,0.25], width = 0.2)
overlaid_ax.set_xlabel(r'$\phi$')
overlaid_ax.set_ylabel(r'$\omega$')
```
![](docs/source/images/line_plot.png?raw=true)


### Text path effects
When embedding the 2D colormap into the plot itself in full or just a part of it, it can happen that the labels and or ticks of the colormap are not visible due to being very similar in color to the picture itself. You can fix this by passing appropriate `path_effects` to add an outline to labels and ticks.

```python
import matplotlib.pyplot as plt
import colorstamps
import matplotlib.patheffects as PathEffects

img = colorstamps.helpers.get_random_data()
rgb, stamp = colorstamps.apply_stamp(
    img[:, :, 0],
    img[:, :, 1],
    "funnel",
    vmin_0=-1.2,
    vmax_0=1.2,
    vmin_1=-1,
    vmax_1=1,
)


fig, axes = plt.subplots(2)

# plot the original image
ax = axes[0]
ax.imshow(rgb)

# show colormap as overlay
overlaid_ax = stamp.overlay_ax(
    ax, lower_left_corner=[0.7, 0.85], width=0.2, path_effects=None
)
overlaid_ax.set_ylabel(r"$\phi$")
overlaid_ax.set_xlabel(r"$\omega$")


# plot the fix
ax = axes[1]
ax.imshow(rgb)

# add path effects to make text more readable
path_effects = [PathEffects.withStroke(linewidth=3, foreground="w")]

# show colormap as overlay
overlaid_ax = stamp.overlay_ax(
    ax, lower_left_corner=[0.7, 0.85], width=0.2, path_effects=path_effects
)
overlaid_ax.set_ylabel(r"$\phi$")
overlaid_ax.set_xlabel(r"$\omega$")
```

![](docs/source/images/path_effects_example.png?raw=true)


# support
Any contributions that would provide additional colormaps for are welcome, as are contributions to increase the functionality.
