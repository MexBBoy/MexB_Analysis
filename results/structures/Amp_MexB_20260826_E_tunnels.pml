load Amp_MexB_20260826_E_tunnels.pdb, scene
hide everything
show cartoon, polymer
color grey80, polymer
color skyblue, polymer and chain E
select tun, resn TUN
alter tun, vdw=b
rebuild
show spheres, tun
set_color cav1, [0.792,0.059,0.757]
color cav1, resn TUN and chain T
set_color cav2, [0.039,0.616,0.627]
color cav2, resn TUN and chain U
set_color cav3, [0.776,0.545,0.235]
color cav3, resn TUN and chain V
show sticks, not polymer and not resn TUN
util.cbay('not polymer and not resn TUN')
bg_color white
set ray_opaque_background, 1
orient
