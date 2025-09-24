TreeIso (plugin)
============

Individual-tree isolation (treeiso) from terrestrial laser scanning
-----------------------------------------------------------------

Zhouxin Xi, Chris Hopkinson (Department of Geography & Environment, University of Lethbridge, Canada)

TreeIso can be utilized to separate individual trees from terrestrial laser scanning point clouds, assigning tree IDs as supplementary scalar fields.

Please cite the following paper if you find this tool helpful:

Xi, Z.; Hopkinson, C. 3D Graph-Based Individual-Tree Isolation (Treeiso) from Terrestrial Laser Scanning Point Clouds. Remote Sens. 2022, 14, 6116. https://doi.org/10.3390/rs14236116

This tool relies on the cut-pursuit algorithm, please also consider citing:
Landrieu, L.; Obozinski, G. Cut Pursuit: Fast Algorithms to Learn Piecewise Constant Functions on General Weighted Graphs. SIAM J. Imaging Sci. 2017, 10, 1724–1766. hal link
Raguet, H.; Landrieu, L. Cut-pursuit algorithm for regularizing nonsmooth functionals with graph total variation. In Proceedings of the International Conference on Machine Learning, Stockholm, Sweden, 10–15 July 2018; Volume 80, pp. 4247–4256.

The cut-pursuit files are licensed under the **MIT** license - see the [LICENSE-MIT](LICENSE-MIT.txt) file for details.

A Matlab version shared via:
https://github.com/truebelief/artemis_treeiso

The development of the treeiso plugin was inspired by the CSF plugin originating from: 
Zhang W, Qi J, Wan P, Wang H, Xie D, Wang X, Yan G. An Easy-to-Use Airborne LiDAR Data Filtering Method Based on Cloth Simulation. Remote Sensing. 2016; 8(6):501


Command line mode
-----------------
Command line is supported. The plugin consists of three stages of segmentation. In the command line mode, please configure the segmentation parameters for each stage to activate the corresponding stage. If no parameters are provided for a specific stage, the segmentation for that stage will not be executed.

Available options
-----------------
<table>
	<tr>
		<th>Command</th>
		<th>Description</th>
	</tr>
	<tr>
		<td><code>-TREEISO</code></td>
		<td>
			<i>Runs the TREEISO plugin</i>
			<p>Optional settings are:</p>
			<ul>
				<li> -LAMBDA1 [value]: Regularization strength for initial segmentation (default: 1.0)</li>
				<li> -K1 [value]: Nearest neighbors to search for initial segmentation(default: 5, minimum>3)</li>
				<li> -DECIMATE_RESOLUTION1 [value]: Decimated point resolution (in m) for initial segmentation (default: 0.05)</li>
				<li> -LAMBDA2 [value]: Regularization strength for intermediate segmentation (default: 20)</li>
				<li> -K2 [value]: Nearest neighbors to search for intermediate segmentation (default: 20)</li>
				<li> -MAX_GAP [value]: Maximum point gap (in m) for intermediate segmentation (default: 2.0)</li>
				<li> -DECIMATE_RESOLUTION2 [value]: Decimated point resolution (in m) for intermediate segmentation (default: 0.1)</li>
				<li> -RHO [value]: Relative height to length ratio (used to detect non-stems for final segmentation) (default: 0.5)</li>
				<li> -VERTICAL_OVERLAP_WEIGHT [value]: Vertical overlapping ratio weight for final segmentation (default: 0.5)</li>
			</ul>
		</td>
	</tr>
</table>

<!-- demo -->
![treeiso_demo](https://user-images.githubusercontent.com/8785889/236364374-5d9f69e0-0877-43b3-9927-f923d65262c1.gif)
<!-- demo -->

## 3rd Party Installation

### [CloudCompare]

You can check out the official build instructions here:
[CloudCompare Build Instructions]

In the CMake, please enable the PLUGIN_STANDARD_QTREEISO option, and you can find the treeiso added to CC after build.

`CMake` Edits:
```diff
# Windows
mkdir build & cd build
- cmake -DCMAKE_PREFIX_PATH=C:\Qt\5.15.2\msvc2019_64 ..
+ cmake -DCMAKE_PREFIX_PATH=C:\Qt\5.15.2\msvc2019_64 -PLUGIN_STANDARD_QTREEISO ..

# macOs
mkdir build && cd build
- cmake -DCMAKE_PREFIX_PATH=/usr/local/opt/qt@5 ..
+ cmake -DCMAKE_PREFIX_PATH=/usr/local/opt/qt@5 -PLUGIN_STANDARD_QTREEISO ..

# Linux
mkdir build && cd build
cmake ..
```

We also have standalone Python and MATLAB versions, but they don’t interact with CloudCompare.

#### Mac OS (Apple Silicon)

The newer **Treeiso** plugin fixes a few minor bugs and is up to 15× faster, but produces similar results as the older one. 

If you’re on a Mac with CC (~2.13.2), you can continue using that version without big issues.

If you are compiling and running locally, add -DCC_MAC_DEV_PATHS to the CMAKE_CXX_FLAGS in the CMAKE group. 
This will look for the plugins in your build directory rather than the application bundle. If you need the shaders as well,
you will have to create a shaders folder in the build directory and copy the shaders you need into it.

For convenience, we provide a 100% self contained build system based on [pixi]:
1. Install pixi from the official website and in the root directory of CC code repository
2. Launch `pixi run build` and then `pixi run CloudCompare`.
3. OPTIONAL: You can also create a portable (relocatable) `.app` by using `pixi run bundle`. You should find the bundle in `.build/install/CloudCompare` directory.
4. OPTIONAL: I haven’t tested this myself, but it might work. Let me know if you give it a try --> Edit the `pixi.toml`, enabling the **TreeISO** plugin by setting `-DPLUGIN_STANDARD_QTREEISO=ON`.

##### Issues

If you run into any issues compiling on macOS, I recommend contacting the CloudCompare founder [@dgirardeau] via [CloudCompare Issues], who knows everything about the releases and cross-platform details. I also hope a recent pre-built mac version can be available, but currently the new plugin is only in the Windows alpha release (2.14.alpha, dated 02/18/2025).

<!-- refs -->
[@dgirardeau]: https://github.com/dgirardeau
[CloudCompare]: https://www.danielgm.net/cc/
[CloudCompare Issues]: https://github.com/CloudCompare/CloudCompare/issues
[CloudCompare Build Instructions]: https://github.com/CloudCompare/CloudCompare/blob/master/BUILD.md
[pixi]: https://pixi.sh/latest/



