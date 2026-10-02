reinitialize
bg_color white
set ray_opaque_background, on
set orthoscopic, on
set antialias, 2
set ray_shadow, off
set depth_cue, off
set ambient, 0.5
set spec_reflect, 0.15
load 1a22.cif, complex
remove hydro
hide everything
show cartoon, complex and polymer.protein
color marine, chain A
color gray75, chain B
set cartoon_fancy_helices, on
set cartoon_transparency, 0.18
select site, chain A and resi 58
show sticks, site
show spheres, site and not name N+C+O
set stick_radius, 0.28, site
set sphere_scale, 0.27, site
color orange, site
orient complex
turn z, -20
zoom complex and polymer.protein, 1.5
set ray_trace_mode, 1
png complex.png, width=1400, height=1100, dpi=450, ray=1
hide cartoon
select local, byres (complex within 7.0 of site)
show sticks, local and not name N+C+O+CA
set stick_radius, 0.10, local and not site
color gray75, local and chain B
color lightblue, local and chain A and not site
show lines, local and name N+C+O+CA
orient local
turn z, 15
zoom site, 6
png site.png, width=1400, height=1000, dpi=450, ray=1
quit
