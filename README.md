<h2>Description of Code</h2>

<h3>scRNA-seq Data</h3>
<ul>
  <li>
    <b>Sample-level processing & QC:</b>
    <code>IRD/scRNA/sample_level_processing/</code>
  </li>
  <li>
    <b>Object creation & merged analyses:</b>
    <code>IRD/scRNA/merged_analysis/</code>
  </li>
  <li>
    <b>Lineage-specific DEG analyses:</b>
    <code>IRD/scRNA/lineage_level_analysis/</code>
  </li>
</ul>

<h3>Xenium Data</h3>
<ul>
  <li>
    <b>Sample-level processing:</b>
    <code>Xenium_to_RDS.R</code>
  </li>
  <li>
    <b>Merged object creation:</b>
    <code>merge_samples.ipynb</code>
  </li>
  <li>
    <b>Cell typing:</b>
    <ul>
      <li><code>manual_cell_typing.ipynb</code></li>
      <li><code>IRD_sketch_annotation.ipynb</code></li>
      <li><code>integrate_manual_and_clustering_celltyping.ipynb</code></li>
    </ul>
  </li>
  <li>
    <b>General Xenium scripts:</b>
    <code>IRD/Xenium/</code>
  </li>
  <li>
    <b>Cell–cell proximity analyses:</b>
    <code>IRD/Xenium/distance_analysis/</code>
  </li>
  <li>
    <b>Radial neighborhood analyses:</b>
    <code>IRD/Xenium/radial_neighborhoods/</code>
  </li>
  <li>
    <b>BANKSY neighborhood analyses:</b>
    <code>IRD/Xenium/BANKSY/</code>
    <br />
    <i>Note:</i> The BANKSY workflow also utilizes code from
    <a href="https://github.com/jwweii/PyBanksy-Harmony" target="_blank">
      PyBanksy-Harmony
    </a>.
  </li>
</ul>

<h3>CODEX Data</h3>
<ul>
  <li>
    <b>Data processing:</b> Utilizes code from
    <a href="https://github.com/ding-lab/ding-lab-spatial/tree/main/multiplex_imaging" target="_blank">
      ding-lab-spatial
    </a>
  </li>
  <li>
    <b>Segmentation, cell typing, and analyses:</b>
    <code>IRD/CODEX/</code>
  </li>
</ul>

<hr />

<h2>Code Used to Generate Figures</h2>

<ul>
  <li>
    <b>Figure 1c–d:</b>
    scRNA-seq and Xenium merging and cell typing scripts described above.
  </li>

  <li>
    <b>Figure 2a–e:</b>
    scRNA-seq cell type abundance comparisons and correlations:
    <ul>
      <li><code>IRD/scRNA/merged_analysis/celltype_abundance_comparison.ipynb</code></li>
      <li><code>IRD/scRNA/merged_analysis/celltype_abundance_correlation.ipynb</code></li>
    </ul>
  </li>

  <li>
    <b>Figure 2f–g:</b>
    Xenium cell type abundance comparisons:
    <code>IRD/Xenium/compare_celltype_fractions.ipynb</code>
  </li>

  <li>
    <b>Figure 3a:</b>
    Cell–cell proximity analyses:
    <code>IRD/Xenium/distance_analysis/</code>
    <br />
    <b>Figure 3b–e:</b>
    Radial neighborhood analyses:
    <code>IRD/Xenium/radial_neighborhoods/</code>
  </li>

  <li>
    <b>Figure 4a–d:</b>
    B cell–specific scRNA-seq analyses:
    <code>IRD/scRNA/lineage_level_analysis/</code>
    <br />
    <b>Figure 4e–f:</b>
    B cell–specific Xenium analyses:
    <code>IRD/Xenium/Bcell_analysis.ipynb</code>
    <br />
    <b>Figure 4g:</b>
    Longitudinal MMRF scRNA-seq analyses:
    <code>IRD/scRNA/MMRF_longitudinal/</code>
  </li>

  <li>
    <b>Figure 5:</b>
    <ul>
      <li>
        Xenium analyses:
        <code>IRD/Xenium/B_interaction_analysis.ipynb</code>
      </li>
      <li>
        scRNA-seq analyses:
        <code>IRD/scRNA/merged_analysis/APRIL_BAFF_analysis.ipynb</code>,
        <code>IRD/scRNA/merged_analysis/APRIL_BAFF_analysis_R.ipynb</code>
      </li>
    </ul>
  </li>

  <li>
    <b>Figure 6:</b>
    Myeloid and mesenchymal stromal cell analyses:
    <ul>
      <li><code>IRD/Xenium/Mye_analysis.ipynb</code></li>
      <li><code>IRD/Xenium/MSC_analysis.ipynb</code></li>
    </ul>
  </li>

  <li>
    <b>Figure 7:</b>
    T cell–related analyses:
    <ul>
      <li><code>IRD/scRNA/lineage_level_analysis/T_NK_subcluster.ipynb</code></li>
      <li><code>IRD/Xenium/Xenium_Tcells_new_celltypes.ipynb</code></li>
      <li><code>IRD/CODEX/merge_codex_slides.ipynb</code></li>
    </ul>
  </li>
</ul>

<h2>Python environment</h2>

```
Python: 3.13.5 | packaged by conda-forge | (main, Jun 16 2025, 08:27:50) [GCC 13.3.0]
Platform: Linux-3.10.0-1160.119.1.el7.x86_64-x86_64-with-glibc2.17
Packages:
  notebook==7.4.5
  yarl==1.20.1
  session_info==1.0.1
  pyarrow==19.0.1
  scikit-learn==1.7.1
  xarray-datatree==0.0.14
  exceptiongroup==1.3.0
  pyparsing==3.2.3
  Brotli==1.1.0
  httpx==0.28.1
  jupyter_core==5.8.1
  optree==0.17.0
  datashader==0.18.2
  comm==0.2.3
  multipledispatch==0.6.0
  soupsieve==2.7
  six==1.17.0
  isoduration==20.11.0
  nbconvert==7.16.6
  webencodings==0.5.1
  uri-template==1.3.0
  zstandard==0.23.0
  babel==2.17.0
  fastjsonschema==2.21.2
  xyzservices==2025.4.0
  arrow==1.3.0
  shapely==2.1.1
  xarray-dataclasses==1.9.1
  idna==3.10
  stdlib-list==0.11.1
  httpcore==1.0.9
  async-lru==2.0.5
  MarkupSafe==3.0.2
  formulaic==1.2.1
  xarray-spatial==0.4.0
  pip==25.1.1
  pyzmq==27.0.2
  mistune==3.1.3
  zict==3.0.0
  importlib_metadata==8.7.0
  overrides==7.7.0
  python-dateutil==2.9.0.post0
  param==2.2.1
  nbformat==5.10.4
  multidict==6.6.3
  matplotlib-inline==0.1.7
  session-info2==0.2.1
  tqdm==4.67.1
  fqdn==1.5.1
  PySide6==6.9.1
  pyproj==3.7.2
  PySocks==1.7.1
  decorator==5.2.1
  platformdirs==4.4.0
  lazy_loader==0.4
  jupyter_client==8.6.3
  tzdata==2025.2
  lz4==4.4.4
  nbclient==0.10.2
  ipykernel==6.30.1
  webcolors==24.11.1
  natsort==8.4.0
  legacy-api-wrap==1.4.1
  fsspec==2023.6.0
  referencing==0.36.2
  tornado==6.5.2
  docrep==0.3.2
  widgetsnbextension==4.0.14
  geopandas==1.1.1
  yq==3.4.3
  aiohappyeyeballs==2.6.1
  parso==0.8.5
  traitlets==5.14.3
  pathlib==1.0.1
  msgpack==1.1.1
  tinycss2==1.4.0
  frozenlist==1.7.0
  seaborn==0.13.2
  contourpy==1.3.3
  asciitree==0.3.3
  bokeh==3.8.0
  hpack==4.1.0
  ipython_pygments_lexers==1.1.1
  h11==0.16.0
  kiwisolver==1.4.9
  pillow==11.3.0
  fasteners==0.19
  lark==1.2.2
  squidpy==1.5.0
  Pygments==2.19.2
  ptyprocess==0.7.0
  colorcet==3.1.0
  imageio==2.37.0
  jupyterlab_widgets==3.0.15
  ipython==9.4.0
  pybind11==3.0.1
  scikit-image==0.25.2
  tomlkit==0.13.3
  pandas==2.3.2
  argon2-cffi==25.1.0
  array-api-compat==1.12.0
  formulaic-contrasts==1.0.0
  sniffio==1.3.1
  interface-meta==1.3.0
  numpy==2.2.6
  attrs==25.3.0
  gmpy2==2.2.1
  aiohttp==3.12.15
  imagecodecs==2024.12.30
  zarr==2.18.7
  munkres==1.1.4
  dask-expr==1.1.21
  jupyter-lsp==2.2.6
  pybind11-global==3.0.1
  numcodecs==0.15.1
  toml==0.10.2
  cached-property==1.5.2
  dask-image==2024.5.3
  jupyterlab_server==2.27.3
  beautifulsoup4==4.13.5
  branca==0.8.1
  pexpect==4.9.0
  jupyter==1.1.1
  PIMS==0.7
  click==8.2.1
  networkx==3.5
  partd==1.4.2
  argcomplete==3.6.2
  umap-learn==0.5.9.post2
  defusedxml==0.7.1
  scipy==1.16.1
  locket==1.0.0
  jupyter_server_terminals==0.5.3
  rpds-py==0.27.0
  typing_extensions==4.15.0
  leidenalg==0.10.2
  dask==2024.12.1
  spatialdata==0.2.5
  websocket-client==1.8.0
  setuptools==80.9.0
  shiboken6==6.9.1
  folium==0.20.0
  mapclassify==2.10.0
  stack_data==0.6.3
  ome-zarr==0.10.2
  more-itertools==10.8.0
  narwhals==2.3.0
  donfig==0.8.1.post1
  rich==14.1.0
  patsy==1.0.1
  jsonschema==4.25.1
  colorama==0.4.6
  rfc3986-validator==0.1.1
  asttokens==3.0.0
  tblib==3.1.0
  scanpy==1.11.4
  cycler==0.12.1
  anndata==0.12.0
  bleach==6.2.0
  cloudpickle==3.1.1
  xarray==2024.9.0
  wcwidth==0.2.13
  sortedcontainers==2.4.0
  mpmath==1.3.0
  Deprecated==1.2.18
  aiosignal==1.4.0
  Jinja2==3.1.6
  omnipath==1.0.8
  cytoolz==1.0.1
  json5==0.12.1
  jupyter-events==0.12.0
  rfc3339_validator==0.1.4
  rfc3987-syntax==1.1.0
  colorist==1.8.6
  sympy==1.14.0
  statannotations==0.7.2
  xmltodict==0.14.2
  spatial_image==1.1.0
  nest_asyncio==1.6.0
  multiscale_spatial_image==1.0.0
  jupyter_server==2.17.0
  statsmodels==0.14.5
  harmonypy==0.0.10
  texttable==1.7.0
  slicerator==1.1.0
  ipywidgets==8.1.7
  matplotlib-scalebar==0.9.0
  types-python-dateutil==2.9.0.20250822
  propcache==0.3.1
  terminado==0.18.1
  toolz==1.0.0
  Send2Trash==1.8.3
  executing==2.2.0
  distributed==2024.12.1
  jsonpointer==3.0.0
  PyYAML==6.0.2
  zipp==3.23.0
  jedi==0.19.2
  h2==4.2.0
  crc32c==2.7.1
  argon2-cffi-bindings==25.1.0
  pure_eval==0.2.3
  joblib==1.5.2
  certifi==2025.8.3
  pynndescent==0.5.13
  typing_utils==0.1.0
  threadpoolctl==3.6.0
  cffi==1.17.1
  jsonschema-specifications==2025.4.1
  jupyter-console==6.6.3
  torch==2.6.0
  pydeseq2==0.5.2
  tifffile==2025.8.28
  pickleshare==0.7.5
  pooch==1.8.2
  pyct==0.5.0
  pandocfilters==1.5.0
  prometheus_client==0.22.1
  llvmlite==0.44.0
  pycparser==2.22
  hyperframe==6.1.0
  validators==0.35.0
  numba==0.61.2
  urllib3==2.5.0
  notebook_shim==0.2.4
  matplotlib==3.10.5
  jupyterlab_pygments==0.3.0
  fonttools==4.59.2
  tomli==2.2.1
  filelock==3.19.1
  psutil==7.0.0
  python-json-logger==2.0.7
  prompt_toolkit==3.0.51
  xarray-schema==0.0.3
  anyio==4.10.0
  requests==2.32.5
  igraph==0.11.9
  markdown-it-py==4.0.0
  pyogrio==0.11.0
  wrapt==1.17.3
  debugpy==1.8.16
  h5py==3.13.0
  jupyterlab==4.4.6
  inflect==7.5.0
  packaging==25.0
  mdurl==0.1.2
  pytz==2025.2
  charset-normalizer==3.4.3
  typeguard==4.4.4
```
<h2>R environment</h2>

```
R version 4.4.3 (2025-02-28)
Platform: x86_64-conda-linux-gnu
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /diskmnt/Users2/jwang/Tools/mamba/envs/cellchat/lib/libopenblasp-r0.3.30.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/Chicago
tzcode source: system (glibc)

attached base packages:
[1] stats4    grid      stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] ggalluvial_0.12.5           DESeq2_1.46.0              
 [3] SummarizedExperiment_1.36.0 Biobase_2.66.0             
 [5] MatrixGenerics_1.18.0       matrixStats_1.5.0          
 [7] GenomicRanges_1.58.0        GenomeInfoDb_1.42.0        
 [9] IRanges_2.40.0              S4Vectors_0.44.0           
[11] BiocGenerics_0.52.0         scProportionTest_0.0.0.9000
[13] enrichR_3.4                 ggpubr_0.6.1               
[15] ComplexHeatmap_2.25.2       plyr_1.8.9                 
[17] lubridate_1.9.4             forcats_1.0.0              
[19] stringr_1.5.1               purrr_1.1.0                
[21] readr_2.1.5                 tibble_3.3.0               
[23] tidyverse_2.0.0             tidyr_1.3.1                
[25] getopt_1.20.4               cowplot_1.2.0              
[27] Matrix_1.7-3                dplyr_1.1.4                
[29] viridis_0.6.5               viridisLite_0.4.2          
[31] reshape2_1.4.4              ggplot2_3.5.2              
[33] Seurat_5.3.0                SeuratObject_5.1.0         
[35] sp_2.2-0                    RColorBrewer_1.1-3         

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22        splines_4.4.3           later_1.4.3            
  [4] polyclip_1.10-7         fastDummies_1.7.5       lifecycle_1.0.4        
  [7] rstatix_0.7.2           doParallel_1.0.17       globals_0.18.0         
 [10] lattice_0.22-7          MASS_7.3-64             backports_1.5.0        
 [13] magrittr_2.0.3          plotly_4.11.0           httpuv_1.6.16          
 [16] sctransform_0.4.2       spam_2.11-1             spatstat.sparse_3.1-0  
 [19] reticulate_1.43.0       pbapply_1.7-4           abind_1.4-5            
 [22] zlibbioc_1.52.0         Rtsne_0.17              WriteXLS_6.8.0         
 [25] circlize_0.4.16         GenomeInfoDbData_1.2.13 ggrepel_0.9.6          
 [28] irlba_2.3.5.1           listenv_0.9.1           spatstat.utils_3.1-5   
 [31] goftest_1.2-3           RSpectra_0.16-2         spatstat.random_3.4-1  
 [34] fitdistrplus_1.2-4      parallelly_1.45.1       codetools_0.2-20       
 [37] DelayedArray_0.32.0     tidyselect_1.2.1        shape_1.4.6.1          
 [40] UCSC.utils_1.2.0        farver_2.1.2            spatstat.explore_3.5-2 
 [43] jsonlite_2.0.0          GetoptLong_1.0.5        progressr_0.15.1       
 [46] Formula_1.2-5           ggridges_0.5.6          survival_3.8-3         
 [49] iterators_1.0.14        foreach_1.5.2           tools_4.4.3            
 [52] ica_1.0-3               Rcpp_1.1.0              glue_1.8.0             
 [55] gridExtra_2.3           SparseArray_1.6.0       withr_3.0.2            
 [58] fastmap_1.2.0           digest_0.6.37           timechange_0.3.0       
 [61] R6_2.6.1                mime_0.13               colorspace_2.1-1       
 [64] scattermore_1.2         tensor_1.5.1            spatstat.data_3.1-8    
 [67] generics_0.1.4          data.table_1.17.8       httr_1.4.7             
 [70] htmlwidgets_1.6.4       S4Arrays_1.6.0          uwot_0.2.3             
 [73] pkgconfig_2.0.3         gtable_0.3.6            lmtest_0.9-40          
 [76] XVector_0.46.0          htmltools_0.5.8.1       carData_3.0-5          
 [79] dotCall64_1.2           clue_0.3-66             scales_1.4.0           
 [82] png_0.1-8               spatstat.univar_3.1-4   tzdb_0.5.0             
 [85] rjson_0.2.23            nlme_3.1-168            curl_7.0.0             
 [88] zoo_1.8-14              GlobalOptions_0.1.2     KernSmooth_2.23-26     
 [91] parallel_4.4.3          miniUI_0.1.2            pillar_1.11.0          
 [94] vctrs_0.6.5             RANN_2.6.2              promises_1.3.3         
 [97] car_3.1-3               xtable_1.8-4            cluster_2.1.8.1        
[100] locfit_1.5-9.12         cli_3.6.5               compiler_4.4.3         
[103] rlang_1.1.6             crayon_1.5.3            future.apply_1.20.0    
[106] ggsignif_0.6.4          stringi_1.8.7           deldir_2.0-4           
[109] BiocParallel_1.40.2     lazyeval_0.2.2          spatstat.geom_3.5-0    
[112] RcppHNSW_0.6.0          hms_1.1.3               patchwork_1.3.2        
[115] future_1.67.0           shiny_1.11.1            ROCR_1.0-11            
[118] igraph_2.1.4            broom_1.0.9
```
