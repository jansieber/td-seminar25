;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "manual"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("scrartcl" "10pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("caption" "labelfont={small}" "textfont={small}") ("subcaption" "") ("tocloft" "") ("amssymb" "") ("graphicx" "") ("amsmath" "") ("helvet" "scaled=0.9") ("beramono" "scaled=0.9") ("mathpazo" "") ("eulervm" "") ("paralist" "") ("typearea" "") ("upquote" "") ("listings" "writefile") ("xcolor" "dvipsnames") ("microtype" "") ("url" "") ("hyperref" "pdftex") ("cleveref" "")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "version"
    "../tools/license"
    "scrartcl"
    "scrartcl10"
    "fontenc"
    "caption"
    "subcaption"
    "tocloft"
    "amssymb"
    "graphicx"
    "amsmath"
    "helvet"
    "beramono"
    "mathpazo"
    "eulervm"
    "paralist"
    "typearea"
    "upquote"
    "listings"
    "xcolor"
    "microtype"
    "url"
    "hyperref"
    "cleveref")
   (TeX-add-symbols
    '("blist" 1)
    '("mlvar" 1)
    '("define" 1)
    '("parm" 1)
    '("file" 1)
    "DDEBIFCODE"
    "ddebifweb"
    "ddebifarx"
    "ddebifwebold"
    "demobase"
    "T"
    "RR"
    "NN"
    "ZZ"
    "CC"
    "defeq"
    "numberthis")
   (LaTeX-add-labels
    "sec:app:get"
    "sec:changes"
    "sec:v32to4a"
    "sec:v32ato32b"
    "sec:v311to32a"
    "sec:v31to311"
    "sec:v3to31"
    "sec:v2to3"
    "int:pocont"
    "sec:intro"
    "sec:minimal"
    "eq:MG"
    "fig:MGbifs"
    "code_struct"
    "struct_pic"
    "explain_dde"
    "sec:ddes"
    "dde"
    "the_dde_type"
    "the_var_equa"
    "A_def"
    "eq:nddepure"
    "sec:dde:stst"
    "eq:deltadef"
    "the_char_eq"
    "sec:dde:psol"
    "sec:dde:hcli"
    "sd_dde"
    "the_dde_type2"
    "the_var_equa2"
    "A_def2"
    "sec:distdelay"
    "DDAE:non-square"
    "dist_delay:disc"
    "dist_delay:simple"
    "chain_delay:yid0"
    "chain_delay:yid"
    "chain_delay:yint"
    "chain_delay:ysum"
    "sec:system:def"
    "defsys:num"
    "defsys:symb"
    "sec:symbolic"
    "sys_def1"
    "example_sys"
    "sec:constrhs"
    "neuron_sys_rhs"
    "neuron_sys_rhs_vec"
    "sec:consttau"
    "sec:symbconst"
    "neuron_symbolic"
    "sec:dirderi"
    "neuron_sys_dirderi"
    "sec:constjac"
    "deri_requested"
    "sys_def2"
    "example_sys2"
    "sec:sdrhs"
    "sec:sdtau"
    "sec:symbsd"
    "sd_symbolic"
    "sec:dirdtau"
    "sec:sdjac"
    "sec:syscond"
    "sec:funcs"
    "sec:symfuncs"
    "sec:dist_delays"
    "sec:distdelay:simple"
    "renewal:demo:math"
    "renewal:dde"
    "renewal:dist"
    "dist_delay:rep"
    "sec:chain:delay"
    "daphnia:r"
    "daphnia:athr"
    "daphnia:bd"
    "daphnia:cd"
    "daphnia:sda"
    "dist:chain:sum"
    "data_structures"
    "sec:funcs:struct"
    "tab:funcs"
    "sec:point:struct"
    "point_structures"
    "sec:stab:struct"
    "stab_structures"
    "sec:method:struct"
    "point_method_structures"
    "meth_stab_struct"
    "continuation_structure"
    "sec:branch:struct"
    "branch_struct"
    "sec:meas:struct"
    "measure_structure"
    "demo"
    "ride-through"
    "sec:demo1:init"
    "sec:demo1:stst"
    "ride1+2_pic"
    "ride3+4_pic"
    "sec:hopf"
    "ride5_pic"
    "ride6_pic"
    "ride7+8_pic"
    "sec:psol"
    "ride9_pic"
    "ride11_pic"
    "ride10+13_pic"
    "ride12_pic"
    "ride14_pic"
    "ride15_pic"
    "demo2"
    "sec:sd:stst"
    "br_stst"
    "sec:sd:hopf"
    "br_hopf"
    "sec:sd:psol"
    "tz_cond"
    "br_ps_sd1"
    "br_ps_sd2"
    "sd_dde_mu"
    "demo3"
    "z"
    "demo3-1"
    "demo3-2"
    "demo3-3"
    "demo3-4+5"
    "demo3-1b"
    "demo3-6"
    "demo3-7"
    "demo3-8"
    "demo3-9"
    "point_manipulation"
    "branch_manipulation"
    "numerical_methods"
    "code_num_methods"
    "determining_systems"
    "determ_stst"
    "determ_fold"
    "determ_hopf"
    "determ_psol"
    "integral_phase_cond"
    "determ_hcli"
    "extra_cond"
    "eq:extra_cond"
    "sys_cond_demo"
    "continuation"
    "root_char_equa_gio_label"
    "lms_method"
    "past_terms"
    "linear_rhs"
    "extract_real_part"
    "extract_imag_part"
    "determ_root"
    "limits_sec"
    "sec:extensions"
    "sec:octave")
   (LaTeX-add-environments
    '("boxit" 1))
   (LaTeX-add-bibliographies)
   (LaTeX-add-saveboxes
    "savepar")
   (LaTeX-add-xcolor-definecolors
    "darkblue"
    "darkred"
    "var"
    "keyword"
    "comment"
    "string"
    "errmsg"))
 :latex)

