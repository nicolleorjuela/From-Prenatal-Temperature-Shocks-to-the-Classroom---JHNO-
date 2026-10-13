*======================================================================
* MASTER.DO - REPLICACION COMPLETA
* Version: octubre 2026
*
* PROYECTO: Efecto de shocks de temperatura en la gestacion y primera
*   infancia (ventanas relativas al mes de nacimiento) sobre puntajes
*   Saber 9 y Saber 11 en Colombia.
*
* QUE HACE ESTE ARCHIVO (en orden):
*   Seccion 1. Construye los paneles (Grado 9 y Grado 11 apilados) para
*              puntaje global, matematicas y lenguaje, con muestra comun.
*   Seccion 2. Construye los pisos termicos a partir de la altitud.
*   Seccion 3. Estima y grafica: niveles Grado 9 vs 11, brechas de genero
*              y urbano-rural, efectos por grupo, y por piso termico.
*   Seccion 4. Construye la base desbalanceada (con y sin Saber 11)
*              uniendo ICFES con los shocks climaticos.
*   Seccion 5. Tobit (robustez por seleccion de quien presenta Saber 11).
*   Seccion 6. Robustez R0-R10.
*   Seccion 7. Mecanismos: NSE, proporcion de ninas, tipo de colegio.
*   Seccion 8. Delta: persistencia formal (efecto 11 menos efecto 9).
*   Seccion 9. Saber 11 multicohorte (2014-II a 2025): escalera de efectos
*              fijos S0-S6, calendario A/B, heterogeneidad y pruebas F1-F3.
*
* INSTRUCCIONES PARA QUIEN LO CORRA:
*  1. Copiar el repositorio.
*  2. Verificar que las bases esten en /data/raw/ :
*       icfes_clima3_19m.dta
*       altitud_robusta.dta
*       altitud_cabecera.csv
*       clima_mensual_1980_2024.dta
*       base_final_sb9_2017_sb11_2018_2024_presentacion.dta
*       icfes_z_2014_2016.dta
*       icfes_z_2017_2019.dta
*       icfes_z_2020_2022.dta
*       icfes_z_2023_2025.dta
*  3. En Stata (version 17 o superior) escribir UNA sola cosa:
*       cd "RUTA/DONDE/QUEDO/LA/CARPETA"
*       do master.do
*     Esa es la unica ruta que se escribe; el resto es automatico.
*  4. Resultados: output/figures (figuras), output/tables (tablas),
*     logs (logs). Los archivos intermedios quedan en data/derived.
*
* Para correr solo una parte, cambiar los interruptores RUN_* de abajo.
*======================================================================
clear all
set more off
version 17
set scheme s1mono          // esquema grafico en blanco y negro
set linesize 120           // ancho de linea en la ventana de resultados

display as result "Replicacion - version octubre 2026 - inicio: `c(current_date)' `c(current_time)'"

*----------------------------------------------------------------------
* 0. RUTAS AUTOMATICAS, CARPETAS Y PAQUETES
*    NO HAY QUE EDITAR NADA AQUI.
*    ROOT = carpeta donde esta Stata (la que contiene master.do).
*    Todas las demas rutas se arman a partir de ROOT.
*----------------------------------------------------------------------
global ROOT "`c(pwd)'"

* Verifica que Stata este parado en la carpeta correcta
capture confirm file "$ROOT/master.do"
if _rc {
    display as error "Haz cd a la carpeta que contiene master.do y vuelve a correr."
    exit 601
}

global RAW "$ROOT/data/raw"            // bases originales (entrada)
global DER "$ROOT/data/derived"        // bases construidas por el codigo
global TMP "$ROOT/data/derived/tmp"    // archivos temporales de coeficientes
global FIG "$ROOT/output/figures"      // figuras PDF
global TAB "$ROOT/output/tables"       // tablas .dta/.csv/.xlsx
global LOG "$ROOT/logs"                // logs de Stata

* Crea las carpetas que no existan
foreach d in data data/raw data/derived data/derived/tmp output output/figures output/tables logs {
    capture mkdir "$ROOT/`d'"
}

* Instala los paquetes necesarios si no estan (reghdfe y ftools)
foreach p in reghdfe ftools {
    capture which `p'
    if _rc ssc install `p', replace
}
capture ftools, compile

*----------------------------------------------------------------------
* INTERRUPTORES Y PARAMETROS
*   RUN_* : 1 = corre esa seccion, 0 = la omite
*----------------------------------------------------------------------
global RUN_FIGS   1     // Secciones 1-3: paneles + figuras (niveles, brechas, grupos, altitud)
global RUN_TOBIT  1     // Seccion 5: Tobit
global RUN_ROBUST 1     // Seccion 6: robustez R0-R10
global RUN_MECH   1     // Seccion 7: mecanismos
global RUN_DELTA  1     // Seccion 8: delta (persistencia)
global RUN_MULTI  1     // Seccion 9: Saber 11 multicohorte
global MULTI_FE   1     // 1 = corre tambien las pruebas de efectos fijos F1-F3 (Seccion 9)

global REHACER 0        // 1 = reconstruye todos los archivos derivados aunque ya existan
global SUF   "_comunFE" // sufijo de los paneles y figuras principales
global TSUF  "_tobit"   // sufijo de las figuras del Tobit
global YB_MIN 2000      // anio de nacimiento minimo (plausible para 9o grado en 2017)
global YB_MAX 2004      // anio de nacimiento maximo
global CORTE 1000       // corte de altitud calor/frio (metros)
global YLIM 0.04        // limite del eje y en figuras del Tobit
global TOBIT_0 1        // Tobit variante 0: solo regresores residualizados
global TOBIT_A 1        // Tobit variante A: + dummies de anio de nacimiento
global TOBIT_B 0        // Tobit variante B: + anio + municipio (LENTO)
global PRUEBA 0         // 1 = usa solo unos departamentos (prueba rapida)
global DEPTOS_PRUEBA "05 08 11 76"   // Antioquia, Atlantico, Bogota, Valle
global NSE_TECNICO 0    // 0 = madre_sup es universitario; 1 = incluye tecnico/tecnologo
global AGE_MIN 15       // edad minima al presentar Saber 11 (multicohorte)
global AGE_MAX 19       // edad maxima al presentar Saber 11 (multicohorte)
global BLOQUES "2014_2016 2017_2019 2020_2022 2023_2025"   // bloques de Saber 11

* Archivos intermedios principales
global FILE_SHOCKS "$DER/clima_mensual_shocks_1980_2024.dta"   // clima con shocks por municipio x mes
global FILE_OUT    "$DER/icfes_clima3_19m_desbal2.dta"         // base desbalanceada con shocks

*----------------------------------------------------------------------
* VENTANAS TEMPORALES (meses relativos al nacimiento)
*   m = antes del nacimiento (m9 = 9 meses antes); 0,1,2,3 = despues.
*   T1: meses -9 a -7 | T2: -6 a -4 | T3: -3 a -1 | T4: 0 a +3 (postnatal)
*   T0 es un placebo de preconcepcion (solo se usa en mecanismos 4.4).
*----------------------------------------------------------------------
global T0 ""
global T1 "m9 m8 m7"
global T2 "m6 m5 m4"
global T3 "m3 m2 m1"
global T4 "0 1 2 3"
global SUF2 ""          // sufijo adicional para nombres de figuras (cambia segun el grupo)

*======================================================================
* LIBRERIA DE PROGRAMAS
* (aqui solo se DEFINEN los programas; no se ejecuta nada todavia)
*======================================================================

*----------------------------------------------------------------------
* need: verifica que una base exista en data/raw; si no, detiene el codigo
*----------------------------------------------------------------------
program define need
    args f
    capture confirm file "$RAW/`f'"
    if _rc {
        display as error "FALTA data/raw/`f'  (ver instrucciones al inicio)"
        exit 601
    }
end

*----------------------------------------------------------------------
* make_bincat: clasifica el shock de temperatura de cada mes en 4 bins
*   1 = shock < -1 DE (frio extremo)
*   2 = [-1, 0) DE   (REFERENCIA)
*   3 = [0, 1) DE
*   4 = >= 1 DE      (calor extremo)
*   Crea bincat_m15 ... bincat_m1, bincat_0 ... bincat_3
*----------------------------------------------------------------------
program define make_bincat
    foreach k in m15 m14 m13 m12 m11 m10 m9 m8 m7 m6 m5 m4 m3 m2 m1 0 1 2 3 {
        gen byte bincat_`k' = .
        replace bincat_`k' = 1 if tshock_`k' <  -1
        replace bincat_`k' = 2 if tshock_`k' >= -1 & tshock_`k' < 0
        replace bincat_`k' = 3 if tshock_`k' >=  0 & tshock_`k' < 1
        replace bincat_`k' = 4 if tshock_`k' >=  1 & !missing(tshock_`k')
    }
end

*----------------------------------------------------------------------
* gen_trim_bins: para una ventana, cuenta cuantos meses cayeron en cada bin
*   args: nombre de la ventana (T1..T4), lista de meses, bin de referencia
*   Crea count_<ventana>_b<bin> (ej. count_T3_b4 = # de meses con calor
*   extremo en la ventana T3). La referencia no se crea (se omite).
*----------------------------------------------------------------------
program define gen_trim_bins
    args trim meses ref
    local bins 1 2 3 4
    local bins_uso : list bins - ref
    foreach b of local bins_uso {
        local suma ""
        foreach k of local meses {
            local suma "`suma' + (bincat_`k'==`b')"
        }
        local suma = subinstr("`suma'"," + ","",1)
        capture drop count_`trim'_b`b'
        quietly gen count_`trim'_b`b' = `suma'
    }
end

*----------------------------------------------------------------------
* get_ctrls_all: arma las listas de controles
*   ctrl_tshock = shocks de los meses que NO estan en ninguna ventana
*                 (entran lineales)
*   ctrl_precip = precipitacion de TODOS los meses (lineal)
*----------------------------------------------------------------------
program define get_ctrls_all
    local todos m15 m14 m13 m12 m11 m10 m9 m8 m7 m6 m5 m4 m3 m2 m1 0 1 2 3
    local enventana $T0 $T1 $T2 $T3 $T4
    local otros : list todos - enventana
    local ctrl_tshock ""
    foreach k of local otros {
        local ctrl_tshock "`ctrl_tshock' tshock_`k'"
    }
    local ctrl_precip ""
    foreach k of local todos {
        local ctrl_precip "`ctrl_precip' precip_`k'"
    }
    c_local ctrl_tshock "`ctrl_tshock'"
    c_local ctrl_precip "`ctrl_precip'"
end

*----------------------------------------------------------------------
* build_globals: reconstruye conteos y globals segun las ventanas actuales
*   VB = variables de interes (conteos por ventana y bin)
*   CT = controles de shocks fuera de ventana
*   CP = controles de precipitacion
*   (se usa en las regresiones de las secciones 5, 6, 7 y 9)
*----------------------------------------------------------------------
program define build_globals
    foreach w in T0 T1 T2 T3 T4 {
        if "${`w'}" != "" {
            gen_trim_bins `w' "${`w'}" 2
        }
    }
    get_ctrls_all
    local vb ""
    foreach w in T0 T1 T2 T3 T4 {
        if "${`w'}" != "" {
            foreach b in 1 3 4 {
                local vb "`vb' count_`w'_b`b'"
            }
        }
    }
    global VB "`vb'"
    global CT "`ctrl_tshock'"
    global CP "`ctrl_precip'"
end

*----------------------------------------------------------------------
* subj: devuelve el nombre de la variable dependiente (outc) segun la
*   asignatura: global -> zpunt_global, math -> zmath, lenguaje -> zlen
*----------------------------------------------------------------------
program define subj
    args s
    if "`s'" == "global" {
        c_local outc "zpunt_global"
    }
    if "`s'" == "math" {
        c_local outc "zmath"
    }
    if "`s'" == "lenguaje" {
        c_local outc "zlen"
    }
end

*----------------------------------------------------------------------
* est_level_all: NIVELES Grado 9 y Grado 11
*   Una sola regresion con las 4 ventanas juntas.
*   Efectos fijos municipio x mes de nacimiento y anio de nacimiento,
*   INTERACTUADOS con grado (once). Errores agrupados por municipio.
*   Guarda, por ventana, el coeficiente y error estandar de cada bin
*   para Grado 9 (coef9, se9) y Grado 11 (coef11, se11).
*   El bin de referencia (2) se guarda con valor 0.
*   Si un coeficiente no se puede estimar, se guarda como missing.
*----------------------------------------------------------------------
program define est_level_all
    args outcome stub
    local wins T1 T2 T3 T4
    local ref 2
    local bins_uso 1 3 4
    foreach w of local wins {
        gen_trim_bins `w' "${`w'}" `ref'
    }
    get_ctrls_all
    local vars_bin ""
    foreach w of local wins {
        foreach b of local bins_uso {
            local vars_bin "`vars_bin' count_`w'_b`b'"
        }
    }
    reghdfe `outcome' c.(`vars_bin' `ctrl_tshock' `ctrl_precip')##i.once, ///
        absorb(muni_month#once year_birth_#once) vce(cluster municipio)
    display as text "N (`outcome', nivel): " e(N)
    foreach w of local wins {
        tempname mh
        postfile `mh' bincode coef9 se9 coef11 se11 using "`stub'_`w'${SUF}.dta", replace
        post `mh' (`ref') (0) (0) (0) (0)
        foreach b of local bins_uso {
            local c9  = .
            local s9  = .
            local c11 = .
            local s11 = .
            * Efecto en Grado 9 = coeficiente base
            capture lincom count_`w'_b`b'
            if _rc == 0 {
                local c9 = r(estimate)
                local s9 = r(se)
            }
            * Efecto en Grado 11 = base + interaccion con once
            capture lincom count_`w'_b`b' + 1.once#c.count_`w'_b`b'
            if _rc == 0 {
                local c11 = r(estimate)
                local s11 = r(se)
            }
            post `mh' (`b') (`c9') (`s9') (`c11') (`s11')
        }
        postclose `mh'
    }
end

*----------------------------------------------------------------------
* est_diff_all: BRECHAS (grupo1 menos grupo0), en Grado 9 y en Grado 11
*   groupvar = mujer (brecha de genero) o urbano (brecha urbano-rural).
*   Una regresion con todo interactuado con grado y con el grupo; los
*   efectos fijos tambien se interactuan con grado y grupo.
*   Se guardan con los mismos nombres (coef9, se9, coef11, se11) para
*   poder usar el mismo programa de graficacion.
*----------------------------------------------------------------------
program define est_diff_all
    args outcome groupvar stub
    local wins T1 T2 T3 T4
    local ref 2
    local bins_uso 1 3 4
    foreach w of local wins {
        gen_trim_bins `w' "${`w'}" `ref'
    }
    get_ctrls_all
    local vars_bin ""
    foreach w of local wins {
        foreach b of local bins_uso {
            local vars_bin "`vars_bin' count_`w'_b`b'"
        }
    }
    reghdfe `outcome' c.(`vars_bin' `ctrl_tshock' `ctrl_precip')##i.once##i.`groupvar', ///
        absorb(muni_month#once#`groupvar' year_birth_#once#`groupvar') vce(cluster municipio)
    display as text "N (`outcome', `groupvar'): " e(N)
    foreach w of local wins {
        tempname mh
        postfile `mh' bincode coef9 se9 coef11 se11 using "`stub'_`w'${SUF}.dta", replace
        post `mh' (`ref') (0) (0) (0) (0)
        foreach b of local bins_uso {
            local d9   = .
            local sd9  = .
            local d11  = .
            local sd11 = .
            * Brecha en Grado 9
            capture lincom c.count_`w'_b`b'#1.`groupvar'
            if _rc == 0 {
                local d9  = r(estimate)
                local sd9 = r(se)
            }
            * Brecha en Grado 11
            capture lincom c.count_`w'_b`b'#1.`groupvar' + c.count_`w'_b`b'#1.once#1.`groupvar'
            if _rc == 0 {
                local d11  = r(estimate)
                local sd11 = r(se)
            }
            post `mh' (`b') (`d9') (`sd9') (`d11') (`sd11')
        }
        postclose `mh'
    }
end

*----------------------------------------------------------------------
* plot_two: figura 2x2 (una celda por ventana) con Grado 9 (azul) y
*   Grado 11 (naranja). Sirve para niveles y para brechas.
*   - Area sombreada: IC al 90% (1.645 errores estandar)
*   - Marcador relleno = significativo al 90%; hueco = no significativo
*   - Eje y fijo en [-0.04, 0.04]; sin titulo general
*   args: archivo de salida, texto del eje y, nota adicional, y los
*   4 archivos de coeficientes (uno por ventana)
*----------------------------------------------------------------------
program define plot_two
    args outfile ytxt note2 f1 f2 f3 f4
    * Une los 4 archivos de coeficientes en uno solo
    tempfile all
    use "`f1'", clear
    gen trim_order = 1
    save `all'
    forvalues t = 2/4 {
        use "`f`t''", clear
        gen trim_order = `t'
        append using `all'
        save `all', replace
    }
    sort trim_order bincode

    * Intervalos de confianza al 90% y marca de significancia
    gen ub9  = coef9  + 1.645*se9
    gen lb9  = coef9  - 1.645*se9
    gen ub11 = coef11 + 1.645*se11
    gen lb11 = coef11 - 1.645*se11
    gen byte sig9  = (lb9  > 0 | ub9  < 0) & !missing(lb9)
    gen byte sig11 = (lb11 > 0 | ub11 < 0) & !missing(lb11)

    local ttitle1 "Trimester 1 (-9 to -7 months)"
    local ttitle2 "Trimester 2 (-6 to -4 months)"
    local ttitle3 "Trimester 3 (-3 to -1 months)"
    local ttitle4 "Postnatal (0 to 3 months)"

    * Un panel por ventana
    forvalues t = 1/4 {
        local yt `"ytitle("")"'
        local leg "legend(off)"
        local xt `"xtitle("")"'
        * Titulo del eje y solo en la columna izquierda
        if `t' == 1 | `t' == 3 {
            local yt `"ytitle("`ytxt'", size(small))"'
        }
        * Titulo del eje x solo en la fila de abajo
        if `t' == 3 | `t' == 4 {
            local xt `"xtitle("Shock Intensity (SD)", size(small))"'
        }
        * Leyenda solo en el primer panel
        if `t' == 1 {
            local leg `"legend(order(3 "Grade 9" 6 "Grade 11") pos(11) ring(0) col(1) size(vsmall) region(lcolor(white) lwidth(vthin)) symxsize(6))"'
        }
        twoway ///
            (rarea ub9 lb9 bincode if trim_order==`t', color("0 114 178%15") lwidth(none)) ///
            (line  coef9 bincode if trim_order==`t', lcolor("0 114 178") lwidth(medthick)) ///
            (scatter coef9 bincode if trim_order==`t' & sig9==0, mcolor("0 114 178") mfcolor(white) mlcolor("0 114 178") msymbol(O) msize(small)) ///
            (scatter coef9 bincode if trim_order==`t' & sig9==1, mcolor("0 114 178") msymbol(O) msize(small)) ///
            (rarea ub11 lb11 bincode if trim_order==`t', color("213 94 0%15") lwidth(none)) ///
            (line  coef11 bincode if trim_order==`t', lcolor("213 94 0") lwidth(medthick)) ///
            (scatter coef11 bincode if trim_order==`t' & sig11==0, mcolor("213 94 0") mfcolor(white) mlcolor("213 94 0") msymbol(D) msize(small)) ///
            (scatter coef11 bincode if trim_order==`t' & sig11==1, mcolor("213 94 0") msymbol(D) msize(small)) ///
            , ///
            title("`ttitle`t''", size(medsmall) color(black)) ///
            xlabel(1 "< -1 SD" 2 "[-1,0)" 3 "[0,1)" 4 ">= 1 SD", labsize(small)) ///
            yscale(range(-0.04 0.04)) ///
            ylabel(-0.04(0.02)0.04, format(%4.2f) labsize(small) angle(0)) ///
            `xt' `yt' ///
            yline(0, lcolor(gs8) lwidth(thin) lpattern(dash)) ///
            `leg' ///
            graphregion(color(white) margin(l+2 r+2 t+1 b+1)) plotregion(color(white) lcolor(gs12)) ///
            name(g`t', replace)
    }

    * Nota al pie (la segunda linea solo si hay nota adicional)
    local n1 "Shaded areas: 90% CI. Reference: [-1,0) SD. Filled markers: significant at 90%."
    local notecmd `"note("`n1'", size(vsmall))"'
    if `"`note2'"' != "" {
        local notecmd `"note("`n1'" "`note2'", size(vsmall))"'
    }
    * Combina los 4 paneles y exporta a PDF
    graph combine g1 g2 g3 g4, rows(2) cols(2) `notecmd' ///
        graphregion(color(white) margin(l+3 r+3)) imargin(medium) xsize(12) ysize(8.5)
    graph export "`outfile'", replace
end

*----------------------------------------------------------------------
* load_panel: carga el panel de una asignatura.
*   Si alt = 1, tambien pega la altitud de la cabecera municipal y crea
*   grupo2 (1 = frio, altitud >= CORTE; 0 = calor).
*----------------------------------------------------------------------
program define load_panel
    args s alt
    use "$DER/panel_`s'_bins4${SUF}.dta", clear
    if `alt' == 1 {
        gen long mun_num = real(municipio)
        merge m:1 mun_num using "$DER/altitud_pisos.dta", keep(match) nogen
        gen byte grupo2 = (alt_cab >= $CORTE) if !missing(alt_cab)
    }
end

*----------------------------------------------------------------------
* run_all: para un subconjunto de la muestra y las 3 asignaturas hace:
*   (a) niveles Grado 9 vs 11, (b) brecha de genero, (c) brecha urbano-rural
*   args: etiqueta (para nombres de archivo), condicion de muestra,
*         texto de la muestra (para la nota), alt (1 = pega altitud)
*----------------------------------------------------------------------
program define run_all
    args tag cond gtxt alt
    local tg ""
    if "`tag'" != "" {
        local tg "_`tag'${SUF2}"
    }
    local ntxt ""
    local stxt ""
    if `"`gtxt'"' != "" {
        local ntxt "Sample: `gtxt'"
        local stxt " Sample: `gtxt'."
    }
    foreach s in global math lenguaje {
        subj `s'

        * (a) Niveles
        load_panel `s' `alt'
        keep if `cond'
        est_level_all `outc' "$TMP/L_`s'`tg'"
        plot_two "$FIG/FIG_gradelevel_`s'`tg'${SUF}.pdf" "Effect (SD)" "`ntxt'" ///
            "$TMP/L_`s'`tg'_T1${SUF}.dta" "$TMP/L_`s'`tg'_T2${SUF}.dta" ///
            "$TMP/L_`s'`tg'_T3${SUF}.dta" "$TMP/L_`s'`tg'_T4${SUF}.dta"

        * (b) y (c) Brechas: genero (mujer) y urbano-rural (urbano)
        foreach g in mujer urbano {
            if "`g'" == "mujer" {
                local gtag "gender"
                local gdef "Gender gap = girls minus boys."
            }
            if "`g'" == "urbano" {
                local gtag "rural"
                local gdef "Urban-rural gap = urban minus rural."
            }
            load_panel `s' `alt'
            keep if `cond'
            est_diff_all `outc' `g' "$TMP/`gtag'_`s'`tg'"
            plot_two "$FIG/FIG_`gtag'_`s'`tg'${SUF}.pdf" "Gap (SD)" "`gdef'`stxt'" ///
                "$TMP/`gtag'_`s'`tg'_T1${SUF}.dta" "$TMP/`gtag'_`s'`tg'_T2${SUF}.dta" ///
                "$TMP/`gtag'_`s'`tg'_T3${SUF}.dta" "$TMP/`gtag'_`s'`tg'_T4${SUF}.dta"
        }
    }
end

*----------------------------------------------------------------------
* run_groups: efectos en NIVELES por grupo, estimados por separado:
*   ninas, ninos, colegios urbanos y colegios rurales.
*   Se omite un grupo si tiene menos de 1,000 observaciones.
*----------------------------------------------------------------------
program define run_groups
    args tag cond gtxt alt
    * Etiquetas de archivo y de nota para cada grupo (0 y 1)
    local lab_mujer0  "ninos"
    local lab_mujer1  "ninas"
    local lab_urbano0 "rural"
    local lab_urbano1 "urbano"
    local txt_mujer0  "Boys"
    local txt_mujer1  "Girls"
    local txt_urbano0 "Rural schools"
    local txt_urbano1 "Urban schools"

    foreach s in global math lenguaje {
        subj `s'
        foreach g in mujer urbano {
            forvalues v = 0/1 {
                local lb "`lab_`g'`v''"
                local full "`lb'"
                if "`tag'" != "" {
                    local full "`tag'_`lb'"
                }
                local full "`full'${SUF2}"
                local gtxt2 "`txt_`g'`v''"
                if `"`gtxt'"' != "" {
                    local gtxt2 "`gtxt' - `txt_`g'`v''"
                }
                load_panel `s' `alt'
                keep if (`cond') & `g' == `v'
                quietly count
                if r(N) < 1000 {
                    display as error "`s' | `full': menos de 1000 obs, se omite"
                    continue
                }
                est_level_all `outc' "$TMP/G_`s'_`full'"
                plot_two "$FIG/FIG_gradelevel_`s'_`full'${SUF}.pdf" "Effect (SD)" "Sample: `gtxt2'" ///
                    "$TMP/G_`s'_`full'_T1${SUF}.dta" "$TMP/G_`s'_`full'_T2${SUF}.dta" ///
                    "$TMP/G_`s'_`full'_T3${SUF}.dta" "$TMP/G_`s'_`full'_T4${SUF}.dta"
            }
        }
    }
end

*----------------------------------------------------------------------
* prep_desbal: prepara la BASE DESBALANCEADA (incluye a quienes
*   presentaron Saber 9 pero no Saber 11) para Tobit, robustez y
*   mecanismos. Crea:
*     tiene11     = 1 si presento Saber 11
*     z9, z11     = puntajes globales estandarizados
*     arriba_50   = 1 si supero la mediana de z11 (no presentados = 0)
*     bincat_*    = bins de shocks
*     muni_month  = grupo municipio x mes de nacimiento
*     depto       = codigo de departamento
*----------------------------------------------------------------------
program define prep_desbal
    use "$FILE_OUT", clear
    foreach v in matemáticas9 lenguaje9 matemáticas11 lenguaje11 {
        capture confirm string variable `v'
        if !_rc {
            destring `v', replace force
        }
    }
    * Cohorte: presentaron Saber 9 (matematicas y lenguaje)
    keep if !missing(matemáticas9, lenguaje9)
    * Anios de nacimiento plausibles
    keep if inrange(year_birth_, $YB_MIN, $YB_MAX)
    * Modo prueba: solo unos departamentos
    if $PRUEBA == 1 {
        gen str2 dpto_str = substr(municipio, 1, 2)
        gen byte en_prueba = 0
        foreach d of global DEPTOS_PRUEBA {
            replace en_prueba = 1 if dpto_str == "`d'"
        }
        keep if en_prueba == 1
        drop dpto_str en_prueba
        display as error "MODO PRUEBA: solo departamentos $DEPTOS_PRUEBA"
    }
    * Llego a Saber 11 (ambos puntajes presentes)
    gen byte tiene11 = !missing(matemáticas11, lenguaje11)
    * Puntaje global = matematicas + lenguaje
    gen punt_global9  = matemáticas9 + lenguaje9
    gen punt_global11 = matemáticas11 + lenguaje11 if tiene11
    * Estandarizar: 9o con toda la cohorte; 11o con quienes presentan
    quietly summarize punt_global9
    gen z9 = (punt_global9 - r(mean)) / r(sd)
    quietly summarize punt_global11 if tiene11
    gen z11 = (punt_global11 - r(mean)) / r(sd) if tiene11
    * Indicador de superar la mediana de 11o (no presentados = 0)
    quietly _pctile z11 if tiene11, p(50)
    gen byte arriba_50 = (tiene11 == 1 & z11 > r(r1))
    make_bincat
    egen muni_month = group(municipio month_birth_)
    gen long mun_num = real(municipio)
    gen depto = real(substr(municipio, 1, 2))
end

*----------------------------------------------------------------------
* post_coefs: guarda coeficientes de la ULTIMA regresion (variables
*   count_*) en un archivo por ventana. Para las figuras de una serie.
*   Opcional: scale multiplica los coeficientes.
*----------------------------------------------------------------------
program define post_coefs
    args stub scale
    if "`scale'" == "" local scale = 1
    foreach w in T1 T2 T3 T4 {
        tempname mh
        postfile `mh' bincode coef se using "`stub'_`w'.dta", replace
        post `mh' (2) (0) (0)
        foreach b in 1 3 4 {
            post `mh' (`b') (_b[count_`w'_b`b']*`scale') (_se[count_`w'_b`b']*`scale')
        }
        postclose `mh'
    }
end

*----------------------------------------------------------------------
* post_coefs_r: igual que post_coefs pero para regresores
*   residualizados (r_count_*), que es lo que usa el Tobit
*----------------------------------------------------------------------
program define post_coefs_r
    args stub scale
    if "`scale'" == "" local scale = 1
    foreach w in T1 T2 T3 T4 {
        tempname mh
        postfile `mh' bincode coef se using "`stub'_`w'.dta", replace
        post `mh' (2) (0) (0)
        foreach b in 1 3 4 {
            post `mh' (`b') (_b[r_count_`w'_b`b']*`scale') (_se[r_count_`w'_b`b']*`scale')
        }
        postclose `mh'
    }
end

*----------------------------------------------------------------------
* run_spec: regresion lineal con efectos fijos exactos
*   (municipio x mes de nacimiento + anio de nacimiento) y guarda
*   los coeficientes para graficar.
*----------------------------------------------------------------------
program define run_spec
    args depvar tag cond
    reghdfe `depvar' $VB $CT $CP if `cond', absorb(muni_month year_birth_) vce(cluster municipio)
    display as text "N (`tag'): " e(N)
    post_coefs "$TMP/tmp_`tag'"
end

*----------------------------------------------------------------------
* run_tobit_fwl: Tobit sobre regresores residualizados (Frisch-Waugh).
*   Guarda dos versiones:
*     _latent: efecto sobre la variable latente (no comparable con OLS)
*     _obs   : efecto marginal sobre el outcome observado (censurado),
*              = efecto latente x P(y* > c)
*----------------------------------------------------------------------
program define run_tobit_fwl
    args tag extra
    tobit y_cens25 $RX `extra' if en_muestra, ll($CLL) vce(cluster municipio)
    display as text "N (`tag'): " e(N) "   censurados Tobit: " e(N_lc)
    post_coefs_r "$TMP/tmp_`tag'_latent"
    * sigma: segun la version de Stata aparece como /sigma o var(e.depvar)
    tempname sg
    capture scalar `sg' = _b[/sigma]
    if _rc {
        scalar `sg' = sqrt(_b[/var(e.y_cens25)])
    }
    * Factor de escala para el efecto marginal observado
    tempvar xb ph
    predict double `xb' if e(sample), xb
    gen double `ph' = normal((`xb' - $CLL)/scalar(`sg')) if e(sample)
    quietly summarize `ph'
    local sc = r(mean)
    display as text "Factor de escala: " %6.3f `sc'
    post_coefs_r "$TMP/tmp_`tag'_obs" `sc'
end

*----------------------------------------------------------------------
* plot_one: figura 2x2 con UNA sola serie (verde) e IC al 90%.
*   Se usa para Tobit y seleccion. Si ylim <= 0, el eje y se ajusta
*   automaticamente a los datos.
*----------------------------------------------------------------------
program define plot_one
    args outfile stub ylim ytxt note
    tempfile all
    use "`stub'_T1.dta", clear
    gen trim_order = 1
    save `all'
    forvalues t = 2/4 {
        use "`stub'_T`t'.dta", clear
        gen trim_order = `t'
        append using `all'
        save `all', replace
    }
    use `all', clear
    sort trim_order bincode
    gen ub = coef + 1.645*se
    gen lb = coef - 1.645*se
    gen byte sig = (lb > 0 | ub < 0) & !missing(lb)

    * Eje automatico si ylim <= 0
    if `ylim' <= 0 {
        quietly summarize ub
        local m1 = r(max)
        quietly summarize lb
        local m2 = abs(r(min))
        local ylim = ceil(max(`m1',`m2')*50)/50
    }
    local y1 = -`ylim'
    local y2 = -`ylim'/2
    local y3 =  `ylim'/2

    local ttitle1 "Trimester 1 (-9 to -7 months)"
    local ttitle2 "Trimester 2 (-6 to -4 months)"
    local ttitle3 "Trimester 3 (-3 to -1 months)"
    local ttitle4 "Postnatal (0 to 3 months)"

    forvalues t = 1/4 {
        local yt `"ytitle("")"'
        local xt `"xtitle("")"'
        if `t' == 1 | `t' == 3 {
            local yt `"ytitle("`ytxt'", size(small))"'
        }
        if `t' == 3 | `t' == 4 {
            local xt `"xtitle("Shock Intensity (SD)", size(small))"'
        }
        twoway ///
            (rarea ub lb bincode if trim_order==`t', color("0 158 115%15") lwidth(none)) ///
            (line  coef bincode if trim_order==`t', lcolor("0 158 115") lwidth(medthick)) ///
            (scatter coef bincode if trim_order==`t' & sig==0, mcolor("0 158 115") mfcolor(white) mlcolor("0 158 115") msymbol(O) msize(small)) ///
            (scatter coef bincode if trim_order==`t' & sig==1, mcolor("0 158 115") msymbol(O) msize(small)) ///
            , ///
            title("`ttitle`t''", size(medsmall) color(black)) ///
            xlabel(1 "< -1 SD" 2 "[-1,0)" 3 "[0,1)" 4 ">= 1 SD", labsize(small)) ///
            yscale(range(`y1' `ylim')) ///
            ylabel(`y1' `y2' 0 `y3' `ylim', format(%5.3f) labsize(small) angle(0)) ///
            `xt' `yt' ///
            yline(0, lcolor(gs8) lwidth(thin) lpattern(dash)) legend(off) ///
            graphregion(color(white) margin(l+2 r+2 t+1 b+1)) plotregion(color(white) lcolor(gs12)) ///
            name(g`t', replace)
    }
    graph combine g1 g2 g3 g4, rows(2) cols(2) ///
        note("Shaded areas: 90% CI. Reference: [-1,0) SD. Filled markers: significant at 90%." ///
             "`note'", size(vsmall)) ///
        graphregion(color(white) margin(l+3 r+3)) imargin(medium) xsize(12) ysize(8.5)
    graph export "`outfile'", replace
end

*----------------------------------------------------------------------
* guarda: agrega a la tabla "resumen" (robustez) los coeficientes de
*   todas las ventanas y bins de la ultima regresion, con una etiqueta
*----------------------------------------------------------------------
program define guarda
    args spec
    foreach w in T1 T2 T3 T4 {
        foreach b in 1 3 4 {
            capture post resumen ("`spec'") ("`w'") (`b') (_b[count_`w'_b`b']) (_se[count_`w'_b`b'])
        }
    }
end

*----------------------------------------------------------------------
* reg_fe: regresion base de robustez (FE municipio x mes + anio,
*   cluster municipio) y guarda los coeficientes con la etiqueta $SEC
*----------------------------------------------------------------------
program define reg_fe
    args depvar cond
    reghdfe `depvar' $VB $CT $CP if `cond', absorb(muni_month year_birth_) vce(cluster municipio)
    guarda "${SEC}_`depvar'"
end

*----------------------------------------------------------------------
* show_t3t4: imprime en pantalla los coeficientes de T3 y postnatal de
*   la ultima regresion (lectura rapida)
*----------------------------------------------------------------------
program define show_t3t4
    foreach w in T3 T4 {
        foreach b in 1 3 4 {
            display as text "  `w' bin `b': b = " %7.4f _b[count_`w'_b`b'] "  se = " %7.4f _se[count_`w'_b`b']
        }
    }
end

*----------------------------------------------------------------------
* gap_g: heterogeneidad por un grupo binario g (madre_sup u oficial)
*   Interactua todos los regresores y los FE con el grupo. Reporta el
*   efecto del grupo 0, del grupo 1, la brecha (1 - 0), su error
*   estandar y p-valor. Omite si algun grupo tiene < 1,000 obs.
*----------------------------------------------------------------------
program define gap_g
    args depvar g cond spec
    * Verifica que ambos grupos tengan suficientes observaciones
    forvalues v = 0/1 {
        quietly count if (`cond') & `g' == `v'
        if r(N) < 1000 {
            display as error "`spec': el grupo `g' == `v' tiene solo " r(N) " obs. Se omite."
            exit
        }
    }
    display as result _n "===== `spec': `depvar' por `g' ====="
    reghdfe `depvar' c.($VB $CT $CP)##i.`g' if (`cond') & !missing(`g'), ///
        absorb(muni_month#`g' year_birth_#`g') vce(cluster municipio)
    foreach w in T1 T2 T3 T4 {
        foreach b in 1 3 4 {
            local v count_`w'_b`b'
            local b0  = .
            local b1  = .
            local gap = .
            local sg  = .
            local pg  = .
            capture quietly lincom `v'
            if !_rc {
                local b0 = r(estimate)
            }
            capture quietly lincom `v' + 1.`g'#c.`v'
            if !_rc {
                local b1 = r(estimate)
            }
            capture quietly lincom 1.`g'#c.`v'
            if !_rc {
                local gap = r(estimate)
                local sg  = r(se)
                local pg  = r(p)
            }
            post resg ("`spec'") ("`g'") ("`w'") (`b') (`b0') (`b1') (`gap') (`sg') (`pg')
            if inlist("`w'","T3","T4") & `b' == 4 {
                display as text "  `w' calor>=1DE | g0: " %7.4f `b0' "  g1: " %7.4f `b1' ///
                    "  brecha: " %7.4f `gap' " (se " %6.4f `sg' ", p = " %5.3f `pg' ")"
            }
        }
    }
    * Prueba conjunta de que la brecha es cero en T3 y T4
    foreach w in T3 T4 {
        display as result "  Prueba conjunta de la brecha en `w' (`spec'):"
        capture noisily test 1.`g'#c.count_`w'_b1 1.`g'#c.count_`w'_b3 1.`g'#c.count_`w'_b4
    }
end

*----------------------------------------------------------------------
* prep_bloque: prepara UN bloque de Saber 11 (multicohorte)
*   - Carga solo las variables necesarias (incluye calendario del colegio)
*   - Construye genero, calendario (A=1, B=0), puntaje (matematicas +
*     lectura critica) y lo estandariza DENTRO de cada periodo de examen
*   - Filtra por edad al presentar (AGE_MIN a AGE_MAX)
*   - Crea bins y conteos por ventana
*   - Guarda el bloque en data/derived/s11_listo2_<bloque>.dta
*----------------------------------------------------------------------
program define prep_bloque
    args b
    use periodo year_birth_ month_birth_ municipio estu_genero ///
        cole_area_ubicacion cole_calendario ///
        punt_matematicas punt_lectura_critica tshock_* precip_* ///
        using "$RAW/icfes_z_`b'.dta", clear
    capture destring periodo, replace force
    capture destring punt_matematicas punt_lectura_critica, replace force
    * Desde 2014-II
    drop if periodo < 20142
    drop if missing(punt_matematicas, punt_lectura_critica, year_birth_, municipio)

    * Genero a partir de estu_genero ("F"/"M"); robusto a "Femenino"/"Masculino"
    capture confirm string variable estu_genero
    if _rc {
        tostring estu_genero, replace force
    }
    gen str1 _g = upper(substr(estu_genero, 1, 1))
    gen byte mujer = (_g == "F") if inlist(_g, "F", "M")
    drop _g estu_genero

    * Calendario del colegio: A = 1, B = 0 (missing si no es A ni B)
    capture confirm string variable cole_calendario
    if _rc {
        tostring cole_calendario, replace force
    }
    gen str1 _c = upper(substr(cole_calendario, 1, 1))
    gen byte cal_A = (_c == "A") if inlist(_c, "A", "B")
    drop _c cole_calendario

    * Puntaje y estandarizacion DENTRO de cada periodo, antes de filtrar por edad
    gen punt_ml = punt_matematicas + punt_lectura_critica
    foreach v in punt_ml punt_matematicas punt_lectura_critica {
        bysort periodo: egen double _m = mean(`v')
        bysort periodo: egen double _s = sd(`v')
        gen float z_`v' = (`v' - _m) / _s
        drop _m _s
    }
    * Edad al presentar
    gen int  anio_ex = floor(periodo/10)
    gen byte edad    = anio_ex - year_birth_
    keep if inrange(edad, $AGE_MIN, $AGE_MAX)

    * Bins y conteos por ventana
    make_bincat
    foreach w in T1 T2 T3 T4 {
        gen_trim_bins `w' "${`w'}" 2
    }
    drop bincat_* punt_*
    compress
    save "$DER/s11_listo2_`b'.dta", replace
end

*----------------------------------------------------------------------
* corre_reg: regresion de la multicohorte. Guarda coeficientes en la
*   tabla "res". Si no se puede estimar, avisa y sigue.
*   Uso: corre_reg depvar [if], spec(nombre) [abs(efectos fijos)]
*   Se llama corre_reg y no "corre" para no chocar con -correlate-.
*----------------------------------------------------------------------
program define corre_reg
    syntax varname [if], SPEC(name) [ABS(string)]
    marksample touse
    if "`abs'" == "" {
        local abs "muni_month year_birth_"
    }
    display "-------- Especificacion: `spec' | `varlist' | FE: `abs' --------"
    capture noisily reghdfe `varlist' $VB $CT $CP if `touse', absorb(`abs') vce(cluster municipio)
    if _rc {
        display as error "La regresion `spec' no se pudo estimar (r(" _rc "))."
        exit
    }
    display as result "  >> `spec': N = " %10.0fc e(N) ///
        "  | T3 calor: b = " %7.4f _b[count_T3_b4] " (se " %6.4f _se[count_T3_b4] ")" ///
        "  | Postnatal calor: b = " %7.4f _b[count_T4_b4] " (se " %6.4f _se[count_T4_b4] ")"
    foreach w in T1 T2 T3 T4 {
        foreach b in 1 3 4 {
            post res ("`spec'") ("`w'") (`b') (_b[count_`w'_b`b']) (_se[count_`w'_b`b'])
        }
    }
end

*======================================================================
* SECCION 1. PANELES POR ASIGNATURA (muestra comun en las 3)
*   Parte de icfes_clima3_19m.dta (estudiantes con Saber 9 y Saber 11).
*   - Define genero (mujer) y zona del colegio (urbano/rural)
*   - Define la MUESTRA COMUN: estudiantes con los 6 puntajes completos
*     (global, matematicas y lenguaje; grados 9 y 11), para que las
*     tres asignaturas se comparen sobre los mismos estudiantes
*   - Estandariza cada puntaje (z) dentro de la muestra comun
*   - Apila Grado 9 y Grado 11 (reshape long) y crea once = 1 en Grado 11
*   Salida: data/derived/panel_{global,math,lenguaje}_bins4_comunFE.dta
*   Solo se reconstruye si no existe o si REHACER = 1.
*======================================================================
if ($RUN_FIGS == 1 | $RUN_DELTA == 1) {
    capture confirm file "$DER/panel_lenguaje_bins4${SUF}.dta"
    if (_rc != 0) | ($REHACER == 1) {
        need icfes_clima3_19m.dta
        use "$RAW/icfes_clima3_19m.dta", clear

        * Puntajes numericos; global 9o = matematicas + lenguaje
        rename punt_global punt_global11
        destring matemáticas9 lenguaje9 matemáticas11 lenguaje11 punt_global11, replace
        gen punt_global9 = matemáticas9 + lenguaje9
        destring punt_global9, replace
        gen id_est = _n          // identificador del estudiante

        * Genero y zona del colegio
        drop if estu_genero == ""
        gen mujer = estu_genero == "Mujer"
        drop if cole_area_ubicacion == ""
        gen rural  = cole_area_ubicacion == "RURAL"
        gen urbano = cole_area_ubicacion == "URBANO"

        * Muestra comun: los 6 puntajes completos
        gen byte muestra_comun = !missing(punt_global9, punt_global11, ///
            matemáticas9, matemáticas11, lenguaje9, lenguaje11)
        count if muestra_comun == 1

        * Bins de shocks y base maestra
        make_bincat
        save "$DER/master_bins4${SUF}.dta", replace

        * Un panel por asignatura
        foreach s in global math lenguaje {
            if "`s'" == "global" {
                local v9 "punt_global9"
                local v11 "punt_global11"
                local z "zpunt_global"
            }
            if "`s'" == "math" {
                local v9 "matemáticas9"
                local v11 "matemáticas11"
                local z "zmath"
            }
            if "`s'" == "lenguaje" {
                local v9 "lenguaje9"
                local v11 "lenguaje11"
                local z "zlen"
            }
            use "$DER/master_bins4${SUF}.dta", clear
            keep if muestra_comun == 1
            * Estandarizacion con toda la muestra comun
            quietly summarize `v9'
            gen `z'9 = (`v9' - r(mean)) / r(sd)
            quietly summarize `v11'
            gen `z'11 = (`v11' - r(mean)) / r(sd)
            keep id_est municipio year_birth_ month_birth_ mujer urbano bincat_* tshock_* precip_* `z'9 `z'11
            * Apila: una fila por estudiante y grado
            reshape long `z', i(id_est) j(grado)
            gen byte once = (grado == 11)           // 1 si es Grado 11
            egen muni_month = group(municipio month_birth_)   // municipio x mes de nacimiento
            xtset id_est grado
            save "$DER/panel_`s'_bins4${SUF}.dta", replace
        }
    }
}

*======================================================================
* SECCION 2. ALTITUD: pisos termicos por altitud de la cabecera
*   Une altitud_robusta.dta con la altitud de la cabecera municipal
*   (altitud_cabecera.csv) y clasifica:
*     piso_cab = 1 Calido (< 1,000 m) | 2 Templado (1,000-2,000 m)
*                3 Frio (>= 2,000 m)
*   cerca_cab marca municipios a menos de 150 m de un corte (sensibles).
*   Salida: data/derived/altitud_pisos.dta
*           output/tables/revision_pisos_cab.xlsx (para revisar a mano)
*   El original altitud_robusta.dta NO se modifica.
*======================================================================
if ($RUN_FIGS == 1) {
    need altitud_robusta.dta
    need altitud_cabecera.csv

    * Altitud de la cabecera municipal
    import delimited "$RAW/altitud_cabecera.csv", clear stringcols(1)
    destring alt_cab, replace force
    gen long mun_num = real(municipio)
    keep mun_num alt_cab
    tempfile cab
    save `cab'

    * Se une con la base de altitud robusta
    use "$RAW/altitud_robusta.dta", clear
    merge 1:1 mun_num using `cab', keep(master match) nogen

    * Clasificacion en pisos termicos
    gen piso_cab = .
    replace piso_cab = 1 if alt_cab <  1000
    replace piso_cab = 2 if alt_cab >= 1000 & alt_cab < 2000
    replace piso_cab = 3 if alt_cab >= 2000 & !missing(alt_cab)
    label define piso 1 "Calido" 2 "Templado" 3 "Frio", replace
    label values piso_cab piso
    count if missing(piso_cab)
    tab piso_cab
    capture tab piso_cab piso_p25, missing
    gen byte cerca_cab = (abs(alt_cab-1000) < 150 | abs(alt_cab-2000) < 150) if !missing(alt_cab)
    save "$DER/altitud_pisos.dta", replace

    * Tabla de revision para chequear municipios cerca de los cortes
    capture export excel municipio nombre depto alt_cab alt_p25 alt_mediana piso_cab piso_p25 ///
        piso_mediana cerca_cab using "$TAB/revision_pisos_cab.xlsx", firstrow(variables) replace
}

*======================================================================
* SECCION 3. FIGURAS PRINCIPALES
*   Bloque A: muestra completa (niveles, brechas y grupos)
*   Bloque B: tres pisos termicos (calido, templado, frio)
*   Bloque C: dos grupos de altitud (calor < CORTE, frio >= CORTE)
*   Cada bloque genera, por asignatura (global, math, lenguaje):
*     - niveles Grado 9 vs 11
*     - brecha de genero y brecha urbano-rural
*     - niveles por grupo (ninas, ninos, urbano, rural)
*   Las figuras salen en output/figures.
*======================================================================
if ($RUN_FIGS == 1) {

    * ---- A. Muestra completa ----
    global SUF2 ""
    run_all    ""  "1" ""  0
    run_groups ""  "1" ""  0

    * ---- B. Tres pisos termicos ----
    global SUF2 "_pisos3"
    run_all    "calido_cab"   "piso_cab == 1" "Warm (<1,000 m)"          1
    run_all    "templado_cab" "piso_cab == 2" "Temperate (1,000-2,000 m)" 1
    run_all    "frio_cab"     "piso_cab == 3" "Cold (>2,000 m)"          1
    run_groups "calido_cab"   "piso_cab == 1" "Warm (<1,000 m)"          1
    run_groups "templado_cab" "piso_cab == 2" "Temperate (1,000-2,000 m)" 1
    run_groups "frio_cab"     "piso_cab == 3" "Cold (>2,000 m)"          1

    * ---- C. Dos grupos (calor / frio) con corte $CORTE ----
    global SUF2 "_corte${CORTE}"

    * Diagnostico: % de estudiantes en "frio" y peso de Bogota segun el corte
    use "$DER/panel_global_bins4${SUF}.dta", clear
    keep if grado == 9
    gen long mun_num = real(municipio)
    merge m:1 mun_num using "$DER/altitud_pisos.dta", keep(match) nogen
    gen bogota = (mun_num == 11001)
    quietly count if !missing(alt_cab)
    local ntot = r(N)
    foreach c in 1200 1500 1800 2000 {
        quietly count if alt_cab >= `c' & !missing(alt_cab)
        local ncold = r(N)
        quietly count if alt_cab >= `c' & bogota == 1
        local nbog = r(N)
        display as text "Corte `c' m: " %4.1f 100*`ncold'/`ntot' "% en frio; Bogota = " %4.1f 100*`nbog'/`ncold' "% del frio"
    }

    run_all    "calor_2g" "grupo2 == 0" "Warm (<$CORTE m)"  1
    run_all    "frio_2g"  "grupo2 == 1" "Cold (>=$CORTE m)" 1
    run_groups "calor_2g" "grupo2 == 0" "Warm (<$CORTE m)"  1
    run_groups "frio_2g"  "grupo2 == 1" "Cold (>=$CORTE m)" 1
    global SUF2 ""
}

*======================================================================
* SECCION 4. BASE DESBALANCEADA CON CHOQUES
*   Necesaria para Tobit, robustez y mecanismos.
*   Incluye estudiantes que presentaron Saber 9 aunque no presentaran
*   Saber 11.
*   Pasos:
*     4a. Construye los shocks de temperatura por municipio x mes:
*         shock = (temp del mes - promedio historico del mes en ese
*         municipio) / desviacion estandar historica
*     4b. Lee la base ICFES, calcula fecha de nacimiento y toma el
*         municipio de Saber 9
*     4c. Une, en un loop, el clima de cada mes de -15 a +3 respecto
*         al mes de nacimiento
*   Salida: data/derived/icfes_clima3_19m_desbal2.dta
*======================================================================
if ($RUN_TOBIT == 1 | $RUN_ROBUST == 1 | $RUN_MECH == 1) {
    capture confirm file "$FILE_OUT"
    if (_rc != 0) | ($REHACER == 1) {

        * ---- 4a. Base climatica con shocks ----
        capture confirm file "$FILE_SHOCKS"
        if _rc != 0 {
            need clima_mensual_1980_2024.dta
            use "$RAW/clima_mensual_1980_2024.dta", clear
            capture rename MPIO_CDPMP municipio
            tostring municipio, replace force
            replace municipio = string(real(municipio), "%05.0f")   // codigo de 5 digitos
            capture confirm variable tm
            if _rc != 0 {
                gen tm = ym(year, month)      // mes calendario
            }
            format tm %tm
            confirm variable municipio year month temp_mean precip_sum
            * Promedio y desviacion historicos por municipio y mes calendario
            bys municipio month: egen temp_meanlong = mean(temp_mean)
            bys municipio month: egen temp_sd       = sd(temp_mean)
            gen temp_shock = (temp_mean - temp_meanlong) / temp_sd
            save "$FILE_SHOCKS", replace
        }

        * ---- 4b. Base ICFES: fecha de nacimiento y municipio ----
        need base_final_sb9_2017_sb11_2018_2024_presentacion.dta
        use "$RAW/base_final_sb9_2017_sb11_2018_2024_presentacion.dta", clear
        confirm variable estu_fechanacimiento

        * Variable con el municipio de Saber 9 (debe haber exactamente una)
        capture ds estu_cod_mcpi*9
        if _rc != 0 | "`r(varlist)'" == "" {
            display as error "No encontre estu_cod_mcpi*9 (municipio de Saber 9)."
            exit 111
        }
        local mun9 `r(varlist)'
        if `: word count `mun9'' != 1 {
            display as error "Se esperaba 1 variable estu_cod_mcpi*9; encontradas: `mun9'"
            exit 111
        }

        * Fecha de nacimiento -> anio, mes y mes calendario
        capture drop year_birth_ month_birth_ tm_birth_ mun_exp municipio
        capture confirm string variable estu_fechanacimiento
        if !_rc {
            gen double ftemp = date(estu_fechanacimiento, "DMY")
            drop estu_fechanacimiento
            rename ftemp estu_fechanacimiento
        }
        format estu_fechanacimiento %td
        gen year_birth_  = year(estu_fechanacimiento)
        gen month_birth_ = month(estu_fechanacimiento)
        gen tm_birth_    = ym(year_birth_, month_birth_)
        format tm_birth_ %tm

        * Municipio (codigo de 5 digitos) = el de Saber 9
        capture confirm string variable `mun9'
        if _rc != 0 {
            tostring `mun9', replace force
        }
        gen str5 municipio = ""
        replace municipio = string(real(`mun9'), "%05.0f") if !missing(real(`mun9'))
        drop if municipio == ""

        * ---- 4c. Merge con clima: meses -15 a +3 respecto al nacimiento ----
        *   Nombres: m15 ... m1 (antes), 0 (mes de nacimiento), 1 2 3 (despues)
        tempfile clk
        forvalues r = -15/3 {
            local a = abs(`r')
            if `r' < 0 {
                local k "m`a'"
            }
            if `r' >= 0 {
                local k "`a'"
            }
            gen tm_rel_`k' = tm_birth_ + `r'
            format tm_rel_`k' %tm
            preserve
                use "$FILE_SHOCKS", clear
                keep municipio tm temp_mean temp_meanlong temp_shock precip_sum
                rename tm            tm_rel_`k'
                rename temp_mean     temp_`k'
                rename temp_meanlong tlong_`k'
                rename temp_shock    tshock_`k'
                rename precip_sum    precip_`k'
                save `clk', replace
            restore
            merge m:1 municipio tm_rel_`k' using `clk', keep(match) nogen
        }
        save "$FILE_OUT", replace
    }
}

*======================================================================
* SECCION 5. TOBIT (robustez por seleccion) - solo puntaje global
*   Problema: solo se observa Saber 11 para quienes lo presentan.
*   Solucion de robustez: se define un outcome censurado por abajo en
*   el percentil 25 de z11 (no presentados y quienes sacan menos del
*   p25 toman el valor del p25).
*   El Tobit no puede absorber ~13,000 dummies de municipio x mes, asi
*   que los regresores se RESIDUALIZAN contra esos efectos fijos
*   (Frisch-Waugh) y el Tobit se corre sobre los residuos.
*   Variantes: 0 (solo residuos), A (+ anio de nacimiento),
*              B (+ anio y municipio; lento).
*   Chequeo: el OLS sobre residualizados debe reproducir el OLS con FE
*   exacto (mismos coeficientes y errores estandar).
*======================================================================
if ($RUN_TOBIT == 1) {
    prep_desbal
    count
    display as result "Muestra final: " r(N)

    * ---- Diagnostico de seleccion (ABK): quien llega a Saber 11 ----
    quietly summarize tiene11
    display as result "Tasa de llegada a Saber 11: " %5.1f 100*r(mean) "%"
    xtile q9 = z9, nq(4)
    tab q9 tiene11, row
    tab year_birth_ tiene11, row

    * ---- Outcomes: P(superar percentil 25/75 de z11); no presentados = 0 ----
    foreach c in 25 75 {
        quietly _pctile z11 if tiene11, p(`c')
        gen byte arriba_`c' = (tiene11 == 1 & z11 > r(r1))
    }
    * ---- Outcome censurado en p25 (para el Tobit) ----
    quietly _pctile z11 if tiene11, p(25)
    local cval = r(r1)
    global CVAL = `cval'
    global CLL  = `cval' + 0.000001      // limite con holgura por redondeo
    gen y_cens25 = cond(tiene11 == 1 & z11 > `cval', z11, `cval')
    quietly count if y_cens25 <= $CLL
    display as text "Proporcion censurada en p25: " %5.1f 100*r(N)/_N "%"

    * Regresores (conteos por ventana/bin y controles)
    build_globals

    * ---- G1. Regresiones lineales con FE exacto ----
    run_spec z9  z9_global  1                  // 9o grado, cohorte completa
    run_spec z11 z11_global 1                  // 11o grado, solo quien presenta
    foreach c in 25 50 75 {
        run_spec arriba_`c' arriba`c'_global 1 // P(superar percentil c), cohorte completa
    }
    run_spec arriba_50 arribaP50_global "tiene11 == 1"   // idem, solo quien presenta
    run_spec tiene11  asiste11_global 1        // P(presenta Saber 11)
    run_spec y_cens25 ycens25_global  1        // OLS sobre outcome censurado

    * ---- G2. Tobit sobre regresores residualizados (FWL) ----
    * 1. Muestra exacta de reghdfe (sin singletons)
    quietly reghdfe y_cens25 $VB $CT $CP, absorb(muni_month year_birth_) vce(cluster municipio)
    gen byte en_muestra = e(sample)
    * 2. Residualizar cada regresor contra municipio x mes + anio de nacimiento
    global RX ""
    foreach x in $VB $CT $CP {
        capture drop r_`x'
        quietly reghdfe `x' if en_muestra, absorb(muni_month year_birth_) residuals(r_`x')
        global RX "$RX r_`x'"
    }
    * 3. Chequeo: OLS sobre residualizados = OLS con FE exacto
    regress y_cens25 $RX if en_muestra, vce(cluster municipio)
    post_coefs_r "$TMP/tmp_ycens25_fwl"
    * 4. Tobit (variantes)
    if $TOBIT_0 == 1 {
        run_tobit_fwl tobit25_fwl0 ""
    }
    if $TOBIT_A == 1 {
        run_tobit_fwl tobit25_fwlA "i.year_birth_"
    }
    if $TOBIT_B == 1 {
        run_tobit_fwl tobit25_fwlB "i.year_birth_ i.mun_num"
    }

    * ---- H. Figuras ----
    * Niveles
    plot_one "$FIG/FIG_z9_global${TSUF}.pdf"  "$TMP/tmp_z9_global"  $YLIM "Effect (SD)" ///
        "Overall Score: grade 9, full cohort (no selection)."
    plot_one "$FIG/FIG_z11_global${TSUF}.pdf" "$TMP/tmp_z11_global" $YLIM "Effect (SD)" ///
        "Overall Score: grade 11, only students who take Saber 11."
    * Seleccion (ABK)
    foreach c in 25 50 75 {
        plot_one "$FIG/FIG_arriba`c'_global${TSUF}.pdf" "$TMP/tmp_arriba`c'_global" $YLIM "Effect on probability" ///
            "Overall Score: P(grade 11 score above percentile `c'), full cohort, non-takers = 0."
    }
    plot_one "$FIG/FIG_arribaP50_global${TSUF}.pdf" "$TMP/tmp_arribaP50_global" $YLIM "Effect on probability" ///
        "Overall Score: P(grade 11 score above percentile 50), only students who take Saber 11."
    plot_one "$FIG/FIG_asiste11_global${TSUF}.pdf" "$TMP/tmp_asiste11_global" $YLIM "Effect on probability" ///
        "P(takes Saber 11), full cohort."
    * Robustez por seleccion: OLS exacto, OLS residualizado, Tobit
    plot_one "$FIG/FIG_ycens25_global${TSUF}.pdf" "$TMP/tmp_ycens25_global" $YLIM "Effect (SD)" ///
        "Overall Score: OLS on score censored at p25, municipality x month FE (exact)."
    plot_one "$FIG/FIG_ycens25_fwl${TSUF}.pdf" "$TMP/tmp_ycens25_fwl" $YLIM "Effect (SD)" ///
        "OLS on censored score, regressors residualized on muni x month + birth year FE (check)."
    foreach v in 0 A B {
        if ${TOBIT_`v'} == 1 {
            plot_one "$FIG/FIG_tobit25_fwl`v'_obs${TSUF}.pdf" "$TMP/tmp_tobit25_fwl`v'_obs" $YLIM ///
                "Effect on censored score (SD)" ///
                "Tobit (variant `v'), residualized regressors, marginal effect on observed outcome."
            plot_one "$FIG/FIG_tobit25_fwl`v'_latent${TSUF}.pdf" "$TMP/tmp_tobit25_fwl`v'_latent" 0 ///
                "Latent effect (SD)" ///
                "Tobit (variant `v'), residualized regressors, latent effect; not comparable in scale with OLS."
        }
    }
}

*======================================================================
* SECCION 6. ROBUSTEZ R0-R10
*   Especificacion base: FE municipio x mes + anio de nacimiento,
*   cluster por municipio. Outcomes: z11 (quienes presentan) y
*   arriba_50 (P(superar p50), cohorte completa).
*   Salidas: logs/log_robustez.log y output/tables/resumen_robustez.csv
*======================================================================
if ($RUN_ROBUST == 1) {
    capture log close _all
    log using "$LOG/log_robustez.log", text replace
    * Tabla donde se acumulan los coeficientes de todas las pruebas
    capture postclose resumen
    postfile resumen str40 spec str4 ventana byte bin double b se using "$TAB/resumen_robustez.dta", replace
    global SEC "R0"

    prep_desbal
    global T4 "0 1 2 3"
    build_globals

    * ---- R0. Diagnostico de identificacion ----
    * Cuanta variacion hay dentro de cada celda municipio x mes?
    display as result _n "===== R0. Diagnostico de identificacion ====="
    tab year_birth_
    egen byte tag_my = tag(muni_month year_birth_)
    bysort muni_month: egen n_anios = total(tag_my)
    tab n_anios                       // anios de nacimiento por celda
    gen byte celda_util = n_anios > 1
    sum celda_util                    // % de estudiantes en celdas que aportan variacion
    * Los choques predicen el anio de nacimiento dentro de la celda? (sin FE de anio)
    reghdfe year_birth_ $VB $CT $CP, absorb(muni_month) vce(cluster municipio)
    testparm $VB                      // si rechaza, hay composicion ligada al choque

    * ---- R1. Censo completo vs panel enlazado (9o) ----
    * Cuanta seleccion hay? z9 se compara en toda la cohorte y solo en quien llega a 11o
    display as result _n "===== R1. Censo completo vs panel enlazado (z9) ====="
    global SEC "R1_censo"
    quietly reg_fe z9 1
    estimates store censo
    global SEC "R1_panel"
    quietly reg_fe z9 "tiene11 == 1"
    estimates store panel
    estimates table censo panel, ///
        keep(count_T3_b1 count_T3_b3 count_T3_b4 count_T4_b1 count_T4_b3 count_T4_b4) b(%7.4f) se(%7.4f)

    * ---- R2. Pruebas conjuntas por ventana ----
    * H0: los tres bins de la ventana son conjuntamente cero
    display as result _n "===== R2. Pruebas conjuntas por ventana ====="
    global SEC "R2_base"
    quietly reg_fe z11 "tiene11 == 1"
    foreach w in T1 T2 T3 T4 {
        test count_`w'_b1 count_`w'_b3 count_`w'_b4
    }
    quietly reg_fe arriba_50 1
    foreach w in T1 T2 T3 T4 {
        test count_`w'_b1 count_`w'_b3 count_`w'_b4
    }

    * ---- R3. Colinealidad ----
    * Correlaciones entre conteos y regresiones con cada ventana sola
    display as result _n "===== R3. Colinealidad ====="
    corr count_T3_b1 count_T3_b4 count_T4_b1 count_T4_b4
    foreach w in T3 T4 {
        reghdfe z11 count_`w'_b1 count_`w'_b3 count_`w'_b4 $CP if tiene11 == 1, ///
            absorb(muni_month year_birth_) vce(cluster municipio)
    }

    * ---- R4. Sin el mes de nacimiento (postnatal = meses 1 a 3) ----
    display as result _n "===== R4. Sin el mes de nacimiento ====="
    global SEC "R4_sin_mes0"
    global T4 "1 2 3"
    build_globals
    reg_fe z11 "tiene11 == 1"
    show_t3t4
    reg_fe arriba_50 1
    show_t3t4
    global T4 "0 1 2 3"          // se restaura la ventana original
    build_globals

    * ---- R5. Cluster por departamento ----
    * (pocos clusters: interpretar con cuidado)
    display as result _n "===== R5. Cluster por departamento ====="
    reghdfe z11 $VB $CT $CP if tiene11 == 1, absorb(muni_month year_birth_) vce(cluster depto)
    guarda "R5_clusterdepto_z11"
    reghdfe arriba_50 $VB $CT $CP, absorb(muni_month year_birth_) vce(cluster depto)
    guarda "R5_clusterdepto_arriba_50"

    * ---- R6. Solo anios de nacimiento modales (2001-2003) ----
    display as result _n "===== R6. Solo anios de nacimiento modales ====="
    global SEC "R6_modales"
    reg_fe z11 "tiene11 == 1 & inlist(year_birth_, 2001, 2002, 2003)"
    show_t3t4
    reg_fe arriba_50 "inlist(year_birth_, 2001, 2002, 2003)"
    show_t3t4

    * ---- R7. Anio de nacimiento x departamento ----
    * Controla la composicion de edades por region
    display as result _n "===== R7. Anio de nacimiento x departamento ====="
    reghdfe z11 $VB $CT $CP if tiene11 == 1, absorb(muni_month year_birth_#depto) vce(cluster municipio)
    guarda "R7_anio_x_depto_z11"
    show_t3t4
    reghdfe arriba_50 $VB $CT $CP, absorb(muni_month year_birth_#depto) vce(cluster municipio)
    guarda "R7_anio_x_depto_arriba_50"
    show_t3t4

    * ---- R8. Especificacion continua ----
    * Un coeficiente por ventana (promedio del shock). No se guarda en
    * el resumen: mirar el log.
    display as result _n "===== R8. Especificacion continua ====="
    egen avg_T1 = rowmean(tshock_m9 tshock_m8 tshock_m7)
    egen avg_T2 = rowmean(tshock_m6 tshock_m5 tshock_m4)
    egen avg_T3 = rowmean(tshock_m3 tshock_m2 tshock_m1)
    egen avg_T4 = rowmean(tshock_0 tshock_1 tshock_2 tshock_3)
    reghdfe z11 avg_T1 avg_T2 avg_T3 avg_T4 $CT $CP if tiene11 == 1, absorb(muni_month year_birth_) vce(cluster municipio)
    reghdfe arriba_50 avg_T1 avg_T2 avg_T3 avg_T4 $CT $CP, absorb(muni_month year_birth_) vce(cluster municipio)

    * ---- R9. Heterogeneidad por clima base (terciles) ----
    * tlong_k = temperatura media historica del mes k; el promedio de 12
    * meses da la temperatura media anual del municipio
    display as result _n "===== R9. Terciles de clima base ====="
    egen double tb_i = rowmean(tlong_m12 tlong_m11 tlong_m10 tlong_m9 tlong_m8 tlong_m7 ///
                               tlong_m6 tlong_m5 tlong_m4 tlong_m3 tlong_m2 tlong_m1)
    bysort municipio: egen double tbase = mean(tb_i)
    xtile tq = tbase, nq(3)           // 1 = frio, 3 = calido
    tabstat tbase, by(tq) stat(min max n)
    forvalues q = 1/3 {
        global SEC "R9_tercil`q'"
        display as result "=== Tercil de clima base `q' (1 = frio) ==="
        reg_fe z11 "tiene11 == 1 & tq == `q'"
        show_t3t4
        reg_fe arriba_50 "tq == `q'"
        show_t3t4
    }

    * ---- R10. Efectos por grupo (no solo brechas) ----
    * Regresion separada por grupo = interactuar todo con el grupo
    display as result _n "===== R10. Efectos por grupo ====="
    capture confirm variable estu_genero
    if !_rc {
        gen byte mujer = (estu_genero == "Mujer") if estu_genero != ""
        forvalues v = 0/1 {
            global SEC "R10_mujer`v'"
            display as result "=== mujer == `v' ==="
            reg_fe z11 "tiene11 == 1 & mujer == `v'"
            show_t3t4
            reg_fe arriba_50 "mujer == `v'"
            show_t3t4
        }
    }
    capture confirm variable cole_area_ubicacion
    if !_rc {
        gen byte urbano = (cole_area_ubicacion == "URBANO") if cole_area_ubicacion != ""
        forvalues v = 0/1 {
            global SEC "R10_urbano`v'"
            display as result "=== urbano == `v' (zona del colegio) ==="
            reg_fe z11 "tiene11 == 1 & urbano == `v'"
            show_t3t4
            reg_fe arriba_50 "urbano == `v'"
            show_t3t4
        }
    }

    * ---- Cierre: resumen en una tabla ----
    postclose resumen
    use "$TAB/resumen_robustez.dta", clear
    gen t = b / se
    format b se t %8.4f
    export delimited using "$TAB/resumen_robustez.csv", replace
    * Vista rapida: calor >= 1 DE en T3 y postnatal
    list spec ventana b se t if bin == 4 & inlist(ventana, "T3", "T4"), sepby(spec) noobs
    log close
}

*======================================================================
* SECCION 7. MECANISMOS
*   4.2 NSE (educacion de la madre): el efecto cambia con los recursos?
*   4.6 Tipo de colegio (oficial vs no oficial): heterogeneidad descriptiva
*   4.4 Proporcion de ninas en la cohorte, con placebo de preconcepcion
*   Incluye verificaciones de que los grupos tienen variacion antes de
*   estimar.
*   Salidas: logs/log_mecanismos_v2.log y output/tables/resumen_mecanismos_v2.csv
*======================================================================
if ($RUN_MECH == 1) {
    capture log close _all
    log using "$LOG/log_mecanismos_v2.log", text replace

    prep_desbal
    global T0 ""
    build_globals

    * ---- Variable: ninas ----
    capture confirm variable estu_genero
    if _rc {
        display as error "No existe estu_genero en la base."
        exit 111
    }
    gen byte mujer = regexm(ustrlower(estu_genero), "^(mujer|f)") if estu_genero != ""
    tab estu_genero mujer, missing          // Mujer = 1, Hombre = 0
    quietly summarize mujer
    if r(min) == r(max) {
        display as error "mujer no tiene variacion: revisar estu_genero."
        exit 459
    }

    * ---- Variable NSE: madre con educacion universitaria ----
    local have_nse 0
    capture confirm variable fami_educacionmadre
    if !_rc {
        local have_nse 1
        capture confirm string variable fami_educacionmadre
        if !_rc {
            gen str80 _madre = fami_educacionmadre
        }
        else {
            decode fami_educacionmadre, gen(_madre)
        }
        replace _madre = ustrlower(ustrtrim(_madre))
        tab _madre, missing
        local patron "universitario"
        if $NSE_TECNICO == 1 {
            local patron "universitario|t.cnico"
        }
        * na / no aplica / no sabe / vacio -> missing; el resto 0/1
        gen byte madre_sup = ustrregexm(_madre, "`patron'") ///
            if !inlist(_madre, "na", "no aplica", "no sabe", "")
        tab _madre madre_sup, missing        // VERIFICA: universitario -> 1; primaria/bachillerato -> 0
        quietly summarize madre_sup
        display as text "madre_sup: media = " %5.3f r(mean) " (esperado ~0.18)"
        if r(min) == r(max) | r(mean) < 0.03 | r(mean) > 0.60 {
            display as error "madre_sup mal definida. Revisa las categorias de _madre."
            exit 459
        }
    }
    else {
        display as error "No esta fami_educacionmadre: se omite 4.2."
    }

    * ---- Variable tipo de colegio: oficial ----
    local have_cole 0
    capture confirm variable cole_naturaleza
    if !_rc {
        local have_cole 1
        capture confirm string variable cole_naturaleza
        if !_rc {
            gen str40 _nat = cole_naturaleza
        }
        else {
            decode cole_naturaleza, gen(_nat)
        }
        replace _nat = ustrlower(ustrtrim(_nat))
        tab _nat, missing
        * Igualdad EXACTA: "no oficial" NO debe contar como oficial
        gen byte oficial = (_nat == "oficial") if inlist(_nat, "oficial", "no oficial")
        tab _nat oficial, missing            // VERIFICA: oficial -> 1; no oficial -> 0
        quietly summarize oficial
        display as text "oficial: media = " %5.3f r(mean) " (esperado ~0.79)"
        if r(min) == r(max) {
            display as error "oficial no tiene variacion."
            exit 459
        }
    }
    else {
        display as error "No esta cole_naturaleza: se omite 4.6."
    }

    * Tabla de resultados de heterogeneidad
    capture postclose resg
    postfile resg str40 spec str12 grupo str4 ventana byte bin double(b0 b1 brecha se_brecha p_brecha) ///
        using "$TAB/resumen_mecanismos_v2.dta", replace

    * 4.2 NSE: compensacion predice efecto postnatal distinto con madre_sup = 1
    if `have_nse' {
        gap_g z11       madre_sup "tiene11 == 1" "NSE_z11"
        gap_g arriba_50 madre_sup "1"            "NSE_arriba50"
    }
    * 4.6 Tipo de colegio
    if `have_cole' {
        gap_g z11       oficial "tiene11 == 1" "COLE_z11"
        gap_g arriba_50 oficial "1"            "COLE_arriba50"
    }
    postclose resg

    * 4.4 Proporcion de ninas (modelo de probabilidad lineal).
    * OJO: mezcla el sexo al nacer con quien llega y se queda en el colegio.
    * Para seleccion al nacer, habria que usar nacimientos DANE por sexo.
    display as result _n "===== 4.4 Proporcion de ninas ====="
    reghdfe mujer $VB $CT $CP, absorb(muni_month year_birth_) vce(cluster municipio)
    foreach w in T1 T2 T3 T4 {
        test count_`w'_b1 count_`w'_b3 count_`w'_b4
    }

    * Version con ventana de preconcepcion (placebo T0: meses -12 a -10)
    display as result _n "===== 4.4 Proporcion de ninas con placebo preconcepcion (T0) ====="
    global T0 "m12 m11 m10"
    build_globals
    reghdfe mujer $VB $CT $CP, absorb(muni_month year_birth_) vce(cluster municipio)
    display as result "--- Prueba conjunta del placebo T0 ---"
    test count_T0_b1 count_T0_b3 count_T0_b4
    global T0 ""
    build_globals

    * Cierre: tabla resumen
    use "$TAB/resumen_mecanismos_v2.dta", clear
    gen t_brecha = brecha / se_brecha
    format b0 b1 brecha se_brecha t_brecha p_brecha %8.4f
    export delimited using "$TAB/resumen_mecanismos_v2.csv", replace
    list spec ventana b0 b1 brecha se_brecha p_brecha if bin == 4 & inlist(ventana,"T3","T4"), sepby(spec) noobs
    log close
}

*======================================================================
* SECCION 8. DELTA: PERSISTENCIA FORMAL (efecto 11o - efecto 9o)
*   Panel apilado (9o y 11o) con efectos fijos de ESTUDIANTE y efectos
*   fijos interactuados con grado. El coeficiente de la interaccion
*   (o_count_*) es delta = efecto en 11o menos efecto en 9o.
*   H0 de persistencia total: delta = 0.
*   Salidas: logs/log_delta.log y output/tables/resumen_delta.csv
*======================================================================
if ($RUN_DELTA == 1) {
    capture log close _all
    log using "$LOG/log_delta.log", text replace

    use "$DER/panel_global_bins4${SUF}.dta", clear
    keep if inrange(year_birth_, $YB_MIN, $YB_MAX)    // mismo filtro que el resto del trabajo
    count
    tab grado

    foreach w in T1 T2 T3 T4 {
        gen_trim_bins `w' "${`w'}" 2
    }
    get_ctrls_all
    local vars_bin ""
    foreach w in T1 T2 T3 T4 {
        foreach b in 1 3 4 {
            local vars_bin "`vars_bin' count_`w'_b`b'"
        }
    }
    local xv `vars_bin' `ctrl_tshock' `ctrl_precip'

    * Interacciones de cada regresor con once (efecto 11o menos efecto 9o)
    local xo ""
    foreach v of local xv {
        capture drop o_`v'
        gen double o_`v' = `v' * once
        local xo "`xo' o_`v'"
    }

    * ---- 1. Niveles por grado (se guardan para la tabla final) ----
    *   g = 0 -> 9o ; g = 1 -> 11o
    foreach g in 0 1 {
        display as result _n "===== Niveles, once == `g' ====="
        quietly reghdfe zpunt_global `vars_bin' `ctrl_tshock' `ctrl_precip' if once == `g', ///
            absorb(muni_month year_birth_) vce(cluster municipio)
        foreach w in T1 T2 T3 T4 {
            foreach b in 1 3 4 {
                scalar b`g'_`w'_`b' = _b[count_`w'_b`b']
                scalar s`g'_`w'_`b' = _se[count_`w'_b`b']
            }
        }
    }

    * ---- 2. Regresion de delta: FE de estudiante + FE x grado ----
    display as result _n "===== delta: panel apilado con FE de estudiante ====="
    reghdfe zpunt_global `xo', absorb(id_est muni_month#once year_birth_#once) vce(cluster municipio)

    tempname ph
    postfile `ph' str4 ventana byte bin double(b9 se9 b11 se11 delta se_delta p_delta lo95 hi95) ///
        using "$TAB/resumen_delta.dta", replace
    foreach w in T1 T2 T3 T4 {
        foreach b in 1 3 4 {
            quietly lincom o_count_`w'_b`b', level(95)
            post `ph' ("`w'") (`b') (scalar(b0_`w'_`b')) (scalar(s0_`w'_`b')) ///
                (scalar(b1_`w'_`b')) (scalar(s1_`w'_`b')) ///
                (r(estimate)) (r(se)) (r(p)) (r(lb)) (r(ub))
        }
    }
    postclose `ph'

    * Pruebas conjuntas de persistencia (H0: delta = 0 para los 3 bins de la ventana)
    foreach w in T1 T2 T3 T4 {
        display as result "--- Cambia el efecto entre 9 y 11 en `w'? ---"
        test o_count_`w'_b1 o_count_`w'_b3 o_count_`w'_b4
    }
    display as result "--- Cambia en alguna ventana? (todos los delta de bins) ---"
    testparm o_count_T?_b?

    * ---- 3. Tabla final ----
    use "$TAB/resumen_delta.dta", clear
    format b9 se9 b11 se11 delta se_delta lo95 hi95 p_delta %8.4f
    export delimited using "$TAB/resumen_delta.csv", replace
    list ventana bin b9 b11 delta se_delta lo95 hi95 p_delta, sepby(ventana) noobs
    log close
}

*======================================================================
* SECCION 9. SABER 11 MULTICOHORTE (2014-II a 2025)
*   Un solo grado (Saber 11), cuatro bloques de 3 anios. Municipio = el
*   de presentacion (igual que el panel).
*
*   ESCALERA DE EFECTOS FIJOS (todas con municipio x mes de nacimiento):
*     S0 : + anio de nacimiento
*     S1 : + anio de nacimiento + periodo de examen
*     S2 : + anio de nacimiento x periodo (misma edad al presentar)
*          <- candidata principal
*     S3 : + anio de nacimiento + departamento x periodo
*     S4 : + anio de nacimiento x periodo + departamento x periodo
*     S5 : S2 + mes de nacimiento x calendario (A/B)
*     S6 : S2 solo con calendario A
*   Cada escalon se corre en tres muestras de anio de nacimiento:
*     todos, 1999-2006 (n9906) y 2000-2004 (n0004, los del panel).
*   Se reportan TODAS las especificaciones; no se escoge por el signo.
*
*   Heterogeneidad (H_*): con S2, por materia, al dia (16-17), sin
*     pandemia, genero y zona del colegio.
*   Pruebas de efectos fijos (F1-F3, si MULTI_FE = 1):
*     F1 : S2 + municipio x periodo
*     F2 : edad fija (16, 17, 18), con FE "plain" y con FE "per" (= S2)
*     F3 : edad 17, nacidos 2000-2004 (plain, per, sin 2020-2021)
*   Salidas: logs/log_multicohorte.log y output/tables/resumen_multicohorte.csv
*======================================================================
if ($RUN_MULTI == 1) {

    * ---- Escalera de efectos fijos ----
    global FE_PLAIN "muni_month year_birth_"
    global FE_S0 "muni_month year_birth_"
    global FE_S1 "muni_month year_birth_ periodo"
    global FE_S2 "muni_month year_birth_#periodo"
    global FE_S3 "muni_month year_birth_ depto#periodo"
    global FE_S4 "muni_month year_birth_#periodo depto#periodo"
    global FE_S5 "muni_month year_birth_#periodo mes_cal"
    global FE_MP "muni_month year_birth_#periodo muni_per"

    * ---- Preparar los cuatro bloques (solo si no existen o REHACER = 1) ----
    *   Se usan archivos s11_listo2_* (con calendario); los s11_listo_*
    *   de versiones anteriores NO se reutilizan.
    foreach b of global BLOQUES {
        capture confirm file "$DER/s11_listo2_`b'.dta"
        if (_rc != 0) | ($REHACER == 1) {
            need icfes_z_`b'.dta
            display "Preparando bloque: `b'"
            prep_bloque `b'
        }
        else {
            display "Se reutiliza el bloque ya preparado: `b'"
        }
    }

    * ---- Juntar los bloques ----
    clear
    foreach b of global BLOQUES {
        append using "$DER/s11_listo2_`b'.dta"
    }
    compress

    * ---- Identificadores y variables adicionales ----
    egen muni_month = group(municipio month_birth_)
    egen muni_per   = group(municipio periodo)      // municipio x periodo de examen
    capture confirm string variable municipio
    if _rc == 0 {
        gen depto = real(substr(municipio, 1, 2))
    }
    else {
        gen depto = floor(municipio/1000)
    }
    gen byte urbano = (cole_area_ubicacion == "URBANO") if cole_area_ubicacion != ""
    * Mes de nacimiento x calendario (solo si el calendario es A o B)
    egen mes_cal = group(month_birth_ cal_A) if !missing(cal_A)

    * ---- Revision de la base ----
    describe, short
    tab anio_ex
    tab edad
    tab year_birth_
    tab mujer, missing         // debe ser cercano a 50% mujeres
    tab urbano, missing
    tab cal_A, missing
    tab periodo cal_A          // correspondencia periodo-calendario
    tabstat z_punt_ml, by(periodo) stat(mean sd n)   // z deberia tener media ~0 por periodo
    tab year_birth_ edad if inrange(year_birth_,2000,2004)

    * ---- Regresores ----
    get_ctrls_all
    local vb ""
    foreach w in T1 T2 T3 T4 {
        foreach b in 1 3 4 {
            local vb "`vb' count_`w'_b`b'"
        }
    }
    global VB "`vb'"
    global CT "`ctrl_tshock'"
    global CP "`ctrl_precip'"

    * ---- Log y archivo de resultados ----
    capture log close _all
    log using "$LOG/log_multicohorte.log", text replace
    capture postclose res
    postfile res str40 spec str4 ventana byte bin double b se using "$TAB/resumen_multicohorte.dta", replace

    * ---- 1. Escalera S0-S4 x tres muestras de anio de nacimiento ----
    foreach n in todos n9906 n0004 {
        if "`n'" == "todos" local cnd "!missing(year_birth_)"
        if "`n'" == "n9906" local cnd "inrange(year_birth_,1999,2006)"
        if "`n'" == "n0004" local cnd "inrange(year_birth_,2000,2004)"
        foreach s in S0 S1 S2 S3 S4 {
            display "============ `s' | nacidos: `n' ============"
            corre_reg z_punt_ml if `cnd', spec(`s'_`n') abs(${FE_`s'})
        }
    }

    * ---- 2. Calendario: S5 (S2 + mes x calendario) y S6 (S2 solo calendario A) ----
    foreach n in todos n0004 {
        if "`n'" == "todos" local cnd "!missing(year_birth_)"
        if "`n'" == "n0004" local cnd "inrange(year_birth_,2000,2004)"
        display "============ S5 | nacidos: `n' ============"
        corre_reg z_punt_ml if `cnd' & !missing(cal_A), spec(S5_`n') abs(${FE_S5})
        display "============ S6 | nacidos: `n' ============"
        corre_reg z_punt_ml if `cnd' & cal_A == 1, spec(S6_`n') abs(${FE_S2})
    }

    * ---- 3. Heterogeneidad y variantes (con S2, muestra completa) ----
    corre_reg z_punt_matematicas,     spec(H_matematicas) abs(${FE_S2})
    corre_reg z_punt_lectura_critica, spec(H_lectura)     abs(${FE_S2})
    corre_reg z_punt_ml if inrange(edad,16,17),        spec(H_al_dia)       abs(${FE_S2})
    corre_reg z_punt_ml if !inlist(anio_ex,2020,2021), spec(H_sin_pandemia) abs(${FE_S2})
    corre_reg z_punt_ml if mujer == 0,  spec(H_hombres) abs(${FE_S2})
    corre_reg z_punt_ml if mujer == 1,  spec(H_mujeres) abs(${FE_S2})
    corre_reg z_punt_ml if urbano == 0, spec(H_rural)   abs(${FE_S2})
    corre_reg z_punt_ml if urbano == 1, spec(H_urbano)  abs(${FE_S2})

    * ---- 4. Pruebas de efectos fijos F1-F3 ----
    if $MULTI_FE == 1 {
        * F1: S2 + municipio x periodo
        corre_reg z_punt_ml, spec(F1_todos) abs(${FE_MP})
        corre_reg z_punt_ml if inrange(year_birth_,2000,2004), spec(F1_n0004) abs(${FE_MP})

        * F2: edad fija, todos los anios de nacimiento (plain vs per)
        foreach e in 16 17 18 {
            corre_reg z_punt_ml if edad == `e', spec(F2_e`e'_plain) abs(${FE_PLAIN})
            corre_reg z_punt_ml if edad == `e', spec(F2_e`e'_per)   abs(${FE_S2})
        }

        * F3: edad 17, nacidos 2000-2004 (los del panel)
        corre_reg z_punt_ml if edad == 17 & inrange(year_birth_,2000,2004), ///
            spec(F3_e17_plain) abs(${FE_PLAIN})
        corre_reg z_punt_ml if edad == 17 & inrange(year_birth_,2000,2004), ///
            spec(F3_e17_per) abs(${FE_S2})
        corre_reg z_punt_ml if edad == 17 & inrange(year_birth_,2000,2004) ///
            & !inlist(anio_ex,2020,2021), spec(F3_e17_sinpand) abs(${FE_PLAIN})
    }

    * ---- Cerrar y guardar resultados ----
    postclose res
    use "$TAB/resumen_multicohorte.dta", clear
    gen t = b/se
    format b se t %8.4f
    export delimited using "$TAB/resumen_multicohorte.csv", replace
    * Vista rapida: calor >= 1 DE en T3 y postnatal
    list spec ventana b se t if bin == 4 & inlist(ventana,"T3","T4"), sepby(spec) noobs
    log close
}

display as result "REPLICACION COMPLETA (version octubre 2026). Fin: `c(current_date)' `c(current_time)'"
display as result "Figuras: output/figures | Tablas: output/tables | Logs: logs"
