#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Helper functions for gc_manager.py

@author: Nick_Singstock
"""
import os
import json
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.core import Lattice
import re
import subprocess

opj = os.path.join
ope = os.path.exists
hartree_to_ev = 27.2114
ang_to_bohr = 1.88973

class helper():
    
    def read_inputs(self, folder, file = 'inputs'):
        # returns a dictionary of inputs
        with open(os.path.join(folder, file), 'r') as f:
            input_text = f.read()
        tags = {}
        for line in input_text.split('\n'):
            if line == '':
                continue
            if line[0] == '#':
                continue
            words = line.split(' ')
            if words[0] in tags:
                if type(tags[words[0]]) != list:
                    tags[words[0]] = [tags[words[0]]]
                tags[words[0]].append(' '.join(words[1:]))
            else:
                tags[words[0]] = ' '.join(words[1:])
        return tags
    
    def write_inputs(self, inputs_dic, folder):
        # expects a dictionary of inputs with either strings or lists of strings as vals
        text = ''
        for k,v in inputs_dic.items():
            if type(v) == list:
                for vi in v:
                    text += k + ' ' + vi + '\n'
            else:
                text += k + ' ' + v + '\n'
        with open(os.path.join(folder, 'inputs'), 'w') as f:
            f.write(text)

    def read_optlog(self, folder, file = 'opt.log', verbose = True):
        try:
            with open(os.path.join(folder, file), 'r') as f:
                opt_text = f.read()
        except:
            print('ERROR: Cannot read opt.log file')
            return False
        if len(opt_text) == 0:
            # no data yet
            return False #[{'energy': 'None', 'force': 'None', 'step': 0}]
        steps = {'1': []}
        step = 1
        not_warned = True
        ionic_step = 0
        for i, line in enumerate(opt_text.split('\n')):
            if 'Convergence Step' in line:
                step = int(line.split()[-1])
                if str(step) not in steps:
                    steps[str(step)] = []
                continue
            if ' Step ' in line or '*Force-consistent' in line:
                continue
            if line == '':
                continue
            var = line.split()
            step_dic = {'opt_method': var[0][:-1],
                        'energy': float(var[3][:-1]) / hartree_to_ev,
                        'force': float(var[4]),
                        'step': ionic_step}
            if float(var[4]) > 10 and not_warned and verbose:
                print('**** WARNING: High forces present in convergence. Check CONTCAR. ****')
                not_warned = False
            ionic_step += 1
            steps[str(step)].append(step_dic)
        return steps
    
    def read_Ecomponents(self, folder):
        try:
            with open(os.path.join(folder, 'Ecomponents'), 'r') as f:
                ecomp_text = f.read()
        except:
            print('ERROR: Cannot read Ecomponents file')
            return False
        ecomp = {'units': 'H'}
        for line in ecomp_text.split('\n'):
            if '----' in line or len(line) == 0:
                continue
            ecomp[line.split()[0]] = float(line.split()[-1])
        return ecomp
    
    def read_outfile(self, folder, contcar = 'None'):
        try:
            if ope(opj(folder, 'out')):
                with open(opj(folder, 'out'), 'r', errors='ignore') as f:
                    out_text = f.read()
            elif ope(opj(folder, 'tinyout')):
                with open(opj(folder, 'tinyout'), 'r', errors='ignore') as f:
                    out_text = f.read()
            else:
                assert False
        except:
            print('Error reading out/tinyout file.')
            return False
        site_data = {}
        record_forces, record_ions = False, False
        el_counter = 0
        net_oxidation = 0
        net_mag = 0
        initial_electrons = None
        final_electrons = None
        
        if contcar != 'None':
            ct_els = list(set([site['label'] for site in contcar['sites']]))
            ct_counter = {el: 0 for el in ct_els}
        
        for li,line in enumerate(out_text.split('\n')):
#            if 'nElectrons' in line and initial_electrons is None:
#                initial_electrons = float(line.split()[1])
            if 'FillingsUpdate' in line:
                try:
                    final_electrons = float(line.split()[4])
                except:
                    pass
            
            if 'Ionic positions in cartesian coordinates' in line:
                record_ions = True
                # reset net oxidation states and magnetization for new group
                net_oxidation = 0
                net_mag = 0
                continue
            if record_ions:
                if 'ion' in line:
                    atom, xpos, ypos, zpos = line.split()[1:5]
                    site_data[str(el_counter)] = {'atom': atom}
                    
                    if contcar != 'None':
                        if el_counter >= len(contcar['sites']):
                            print('Error: different number of sites in contcar and out file.')
                            return False
                        try:
                            el_index = [i for i, site in enumerate(contcar['sites']) if atom == site['label']]
                            ct_index = el_index[ct_counter[atom]]
                            ct_site = contcar['sites'][ct_index]
                            ct_counter[atom] += 1
                        except:
                            print('Error reading out file sites.')
                            return False
                        site_data[str(el_counter)]['positions'] = ct_site['xyz'] 
                        site_data[str(el_counter)]['contcar_index'] = ct_index
                        
                        if 'properties' in ct_site and 'selective_dynamics' in ct_site['properties']:
                            site_data[str(el_counter)]['selective_dynamics'] = (
                                        ct_site['properties']['selective_dynamics'])
                    el_counter += 1
                else: #if line == '' or 'ion' not in line:
                    record_ions = False
                    el_counter = 0
                    if contcar != 'None':
                        ct_counter = {el: 0 for el in ct_els}
            
            if 'Forces in Cartesian coordinates' in line:
                record_forces = True
                continue
            if record_forces:
                if 'force' in line:
                    atom, xfor, yfor, zfor = line.split()[1:5]
                    if atom != site_data[str(el_counter)]['atom']:
                        print('Error reading out file.')
                        return False
                    site_data[str(el_counter)]['forces'] =  [float(xfor), float(yfor), float(zfor)]
                    el_counter += 1
                else: #if line == '' or 'force' not in line:
                    record_forces = False
                    el_counter = 0
            
            if 'magnetic-moments' in line:
                el = line.split()[2]
                mags = line.split()[3:]
                mag_counter = 0
                for site,sv in site_data.items():
                    if sv['atom'] == el:
                        try:
                            site_data[site]['mag_mom'] = float(mags[mag_counter])
                            net_mag += float(mags[mag_counter])
                            mag_counter += 1
                        except:
                            print('Error reading magnetic moments')
                            return False
            if 'oxidation-state' in line:
                el = line.split()[2]
                oxis = line.split()[3:]
                oxi_counter = 0
                for site,sv in site_data.items():
                    if sv['atom'] == el:
                        try:
                            site_data[site]['oxi_state'] = float(oxis[oxi_counter])
                            net_oxidation += float(oxis[oxi_counter])
                            oxi_counter += 1
                        except:
                            print('Error reading oxidation states')
                            return False
        return site_data, net_oxidation, net_mag, final_electrons, initial_electrons
    
    def read_out_struct(self, folder, site_data = False):
        try:
            with open(opj(folder, 'out'), 'r', errors='ignore') as f:
                out_text = f.read()
        except:
            assert False, 'Error reading out file.'
                
        record_ions = False
        record_lattice = False
        lattice = np.zeros((3,3))
        no_atoms = True
        no_lattice = True
        if site_data:
            spdata = {'selective_dynamics': []}
        for line in out_text.split('\n'):
            # get lattice
            if '# Lattice vectors:' in line:
                record_lattice = True
                lattice_track = 0
                continue
            if record_lattice:
                if 'unit cell volume' in line:
                    record_lattice = False
                    lattice = lattice.T
                    continue
                if 'R =' in line:
                    continue
                no_lattice = False
                lattice[lattice_track] = np.array([float(x)/ang_to_bohr for x in line.split()[1:-1]])
                lattice_track += 1
            
            # get ions
            if 'Ionic positions in cartesian coordinates' in line:
                record_ions = True
                coords = []
                species = []
                continue
            if record_ions:
                if 'ion' not in line:
                    record_ions = False
                    continue
                species.append(line.split()[1])
                coords.append([float(x)/(ang_to_bohr) for i,x in enumerate(line.split()[2:5])])
                no_atoms = False
                if site_data:
                    spdata['selective_dynamics'] += int(line.split()[-1])
        
        if no_lattice or no_atoms:
            assert False, 'No lattice or no atoms. Out file structure not read correctly.'
        if site_data:
            return Structure(lattice, species, coords, coords_are_cartesian=True, site_properties=spdata)
        return Structure(lattice, species, coords, coords_are_cartesian=True)
                
        

    def read_contcar(self, folder):
        st = Structure.from_file(os.path.join(folder, 'CONTCAR'))
        return st.as_dict()
            
    def get_force(self, steps):
        return steps[-1]['force']
    
    def get_energies(self, steps):
        return [s['energy'] for s in steps]
    
    def check_convergence(self, folder, inputs, steps):
        current_step = list(steps.keys())[-1]
        fmax = 'None'
        econv = 'None'
        max_steps = 'None'
        if ope(opj(folder, 'convergence')):
            conv = self.read_convergence(folder)
            last_step = list(conv.keys())[-1]
            if current_step != last_step: # convergence is on last step
                return False
            if 'fmax' in conv[current_step]:
                fmax = float(conv[current_step]['fmax'])
            if 'econv' in conv[current_step]:
                econv = float(conv[current_step]['econv'])
            if 'max_steps' in conv[current_step]:
                max_steps = int(conv[current_step]['max_steps'])
        
        if fmax == 'None':
            fmax = float(inputs['fmax'])
        if econv == 'None':
            econv = float(inputs['econv']) if 'econv' in inputs else 'None'
        if max_steps == 'None':
            max_steps = int(inputs['max_steps'])
        
        # check if convergence is complete based on force or energy convergence criteria (may change per step)
        force = self.get_force(steps[current_step]) if len(steps[current_step]) > 0 else 1e6
        energies = self.get_energies(steps[current_step]) if len(steps[current_step]) > 0 else []
        nsteps = len(steps[current_step])
        
        if force <= fmax:
            return True # force based convergence
        elif econv != 'None' and len(energies) > 2 and np.abs(energies[-1]-energies[-2]) < econv:
            return True # energy based convergence
        elif nsteps > 0 and max_steps <= 1:
            return True # convergence based on single point calc (no conv for hitting max_steps otherwise)
        return False
    
    def get_bias(self, bias_str):
        if bias_str in ['No_bias','None']:
            return 'No_bias'
        return float(bias_str[:-1])
    
    def get_bias_str(self, bias):
        if bias == 'No_bias':
            return 'No_bias'
        return '%.2f'%bias + 'V'
    
    def read_data(self, folder):
        # currently reads inputs, opt_log for energies, and CONTCAR. Also checks convergence based on forces
        # reads out file for oxidation states and magentic moments
        inputs = self.read_inputs(folder)
        opt_steps = self.read_optlog(folder)
        
        if opt_steps == False:
            return {'opt': 'None', 'inputs': inputs, 'converged': False, 'root': folder,
                    'final_energy': 'None', 'contcar': 'None', 'current_force': 'None'}
        # check if calc has high forces
        current_step = list(opt_steps.keys())[-1]
        current_force = opt_steps[current_step][-1]['force'] if len(opt_steps[current_step]) > 0 else 'None'
        current_energy = opt_steps[current_step][-1]['energy'] if len(opt_steps[current_step]) > 0 else 'None'
        if current_force != 'None' and current_force > 10:
            print('**** WARNING: High forces (> 10) in current step! May be divergent. ****')
        
        # check for convergence
        conv = self.check_convergence(folder, inputs, opt_steps)
        contcar = 'None'
        if 'CONTCAR' in os.listdir(folder):
            contcar = self.read_contcar(folder)
        if not conv:
            return {'opt': opt_steps, 'current_force': current_force, 'current_step': current_step,
                'inputs': inputs, 'current_energy': current_energy,
                'converged': conv, 'contcar': contcar, 
                'final_energy': 'None', 'energy_units': 'H', 'root': folder}
        
        ecomp = self.read_Ecomponents(folder)
        eigStats = self.read_eigStats(folder)
        atom_forces = self.read_forces(folder)
        
        self.make_tinyout(folder)
        
        if ope(opj(folder, 'convergence')):
            convergence = self.read_convergence(folder)
        else:
            convergence = {}
        
        out_sites = self.read_outfile(folder, contcar)
        if out_sites == False:
            sites = {}
            net_oxi = 'None'
            net_mag = 'None'
            nfinal = 'None'
        else:
            sites = out_sites[0]
            net_oxi = out_sites[1]
            net_mag = out_sites[2]
            nfinal = out_sites[3]
        return {'opt': opt_steps, 'current_force': current_force, 'current_step': current_step,
                'inputs': inputs, 'Ecomponents': ecomp, 'current_energy': current_energy,
                'Ecomp_energy': ecomp['F'] if 'F' in ecomp else (ecomp['G'] if 'G' in ecomp else 'None'),
                'converged': conv, 'contcar': contcar, 'nfinal': nfinal,
                'final_energy': 'None' if not conv else current_energy,
                'site_data': sites, 'net_oxidation': net_oxi, 'net_magmom': net_mag,
                'convergence_file': convergence, 'eigStats': eigStats, 'energy_units': 'H',
                'atom_forces': atom_forces, 'root': folder}
        
    def get_neb_data(self, folder, bias):
        # reads neb folder and returns data as a dictionary
        inputs = self.read_inputs(folder)
        opt_steps = self.read_optlog(folder, 'neb.log')
        if opt_steps == False:
            return {'opt': 'None', 'inputs': inputs, 'converged': False,
                    'final_energy': 'None', 'images': {}}
        # check if calc has high forces
        current_step = list(opt_steps.keys())[-1]
        current_force = opt_steps[current_step][-1]['force'] if len(opt_steps[current_step]) > 0 else 'None'
        current_energy = opt_steps[current_step][-1]['energy'] if len(opt_steps[current_step]) > 0 else 'None'
        if current_force != 'None' and current_force > 10:
            print('**** WARNING: High forces (> 10) in current step! May be divergent. ****')
        # check for convergence
        conv = self.check_convergence(folder, inputs, opt_steps)
        images = {}
        for root, folders, files in os.walk(folder):
            # look at subfolders
            if 'CONTCAR' not in files or 'opt.log' not in files:
                continue
            image_num = root.split(os.sep)[-1]
            contcar = self.read_contcar(root)
            ecomp = self.read_Ecomponents(root)
            energy = ecomp['F'] if 'F' in ecomp else (ecomp['G'] if 'G' in ecomp else 'None')
            out_sites = self.read_outfile(root, contcar)
            if out_sites == False:
                sites = {}
                net_oxi = 'None'
                net_mag = 'None'
                nfinal = 'None'
            else:
                sites = out_sites[0]
                net_oxi = out_sites[1]
                net_mag = out_sites[2]
                nfinal = out_sites[3]
            images[image_num] = {'contcar': contcar, 'energy': energy, 'Ecomponents': ecomp, 'nfinal': nfinal,
                                 'site_data': sites, 'net_oxidation': net_oxi, 'net_magmom': net_mag}
        return {'opt': opt_steps, 'current_force': current_force, 'current_step': current_step,
                'inputs': inputs, 'current_energy': current_energy,
                'converged': conv,
                'final_energy': 'None' if not conv else current_energy,
                'images': images, 'path_energy': {k:v['energy'] for k,v in images.items()}}
        
    def get_sites(self, line):
        sites = line.split('[')[-1].replace(']','').split(', ')
        for i,s in enumerate(sites):
            # sites can be a type of site (ie hollow), 'all', 'center', int, element symbol, or (x,y,z) tuple
            try:
                sites[i] = int(s)
                continue
            except:
                pass
            if '(' in s and ')' in s:
                try:
                    xyz = s.replace(')','').replace('(','').split(',')
                    xyz = tuple([float(x) for x in xyz])
                    sites[i] = xyz
                except:
                    pass
            # strings are already treated as themselves in initial split
#            if s in ['All', 'Hollow', 'Ontop', 'Bridge', 'all', 'hollow', 'ontop', 'bridge', 'center']:
#                sites[i] = s
#                continue
        return sites

    def check_surface(self, file, dist = 6):
        # check that surface seems reasonable for manager to handle
        st = Structure.from_file(file)
        # check that z-direction is long distance
        if not (st.lattice.c > st.lattice.a and st.lattice.c > st.lattice.b):
            print('Bad Surface: z-direction is not longest side, skipping: '+file)
            return False
        # check that 10A (dist) of space is available in the vacuum above surface
        coords = st.cart_coords
        max_lens = [np.max(coords[:,0]), np.max(coords[:,1]), np.max(coords[:,2])]
#        if max_lens[2] <= max_lens[0] or max_lens[2] <= max_lens[1]:
#            print('Bad Surface: surface is wider than it is tall, skipping: '+file)
#            return False
        if st.lattice.c - max_lens[2] < dist:
            print('Bad Surface: Needs at least '+str(dist)+'A of vacuum space above surface, skipping: '+file)
            return False
        return True

    def read_bias(self, line):
        bias_str = line.split('[')[-1].replace(']','').split(', ')
        return ['No_bias' if x in ['None','none','No_bias'] else float(x) for x in bias_str]
    
    def update_run_new(self, run_new):
        for root in run_new:
            inputs = self.read_inputs(root)
            inputs['restart'] = 'False'
            self.write_inputs(inputs, root)
    
    def make_tinyout(self, root):
        # only run for newly converged calcs
        # creates 'tinyout' which only has last 'out' calculation section
        if ope(opj(root, 'tinyout')):
            print('tinyout already exists. Skipping.')
            return 
        try:
            with open(opj(root, 'out'), 'r') as f:
                out = f.read()
        except:
            print('ERROR: tinyout/out file not found for converged calc.')
            return 
        tinyout = ''
        for line in out.split('\n'):
            if 'Start date and time' in line:
                tinyout = line + '\n'
            else:
                tinyout += line + '\n'
        with open(opj(root, 'tinyout'), 'w') as f:
            f.write(tinyout)

    def analyze_data(self, all_data):
        '''
        Creates analysis.json
        
        DONE: make analysis.csv (viewable with excel)
        '''
        #print('Data analysis not yet available. Please contact Nick to add.')
        
        # scan through all entries (folders within subdirs)
        self.mols = {}
        bulks = {}
        
#        assert 'O2' in all_data, 'METAERROR: O2 not in all_data'
        
        for entry, entryv in all_data.items():
            if entry == 'converged':
                continue
            # 1. deal with molecules (move to dict for now)
            if 'surf' not in entryv and 'bulk' not in entryv:
                if entry not in self.reference_molecules():
                    print('Entry not recognized: '+entry)
                    continue
                self.mols[entry] = entryv
            # 2. deal with bulk systems (move to dict for now)
            if 'surf' not in entryv and 'bulk' in entryv:
                bulks[entry] = entryv
                
#        print(self.mols.keys())
        
        # create analysis dictionary
        analysis = {}
        for entry, data in all_data.items():
            if entry == 'converged':
                continue
            if 'adsorbed' not in data or 'surf' not in data:
                continue # need surface and adsorbates converged to do analysis on system 
            
            print('Analysis of system: '+entry)
            
            analysis[entry] = {'data': data}
            if entry.split('_')[0] in bulks:
                analysis[entry]['data']['bulk'] = bulks[entry.split('_')[0]]
            else:
                analysis[entry]['data']['bulk'] = 'None'
            
            # run system analysis function to get metrics 
            analysis[entry]['analyzed'] = self.system_analysis(analysis[entry]['data'])
        
        # save data
        file = 'results/analyzed.json'
        with open(file,'w') as f:
            json.dump(analysis, f)
        
        # make analysis.csv
        self.csv_analysis(analysis)
        
        return None
    
    def order_bias(self, bias_list):
        bias_floats = [float(b.replace('V','')) for b in bias_list if b not in ['No_bias']]
        bias_floats.sort()
        bias_sort = [b for bf in bias_floats for b in bias_list if 
                     b != 'No_bias' and float(b.replace('V','')) == bf]
        if 'No_bias' in bias_list:
            bias_sort.append('No_bias')
        assert len(bias_sort) == len(bias_list), 'METAERROR: Sorted bias list is incorrect length.'
        return bias_sort
    
        
    
    def csv_analysis(self, analysis):
        string = 'Ref. Type/Surf,Atom/Mol/Ads,Site_Number,Site_Atom,Ref./Ads. Energies at Biases (eV)\n'
        
        # add ref atoms
        biases = self.order_bias(list(set([ b for mol in self.mols for b in self.mols[mol] ])))
        string += 'Ref. Atoms,,,,' + ','.join(biases) + '\n\n'
        ref_atoms = {}
        for mol, bias_data in self.mols.items():
            for bias in bias_data:
                ra = self.ref_energies(self.formula_to_dic(mol), bias)
                if ra is None:
                    continue
                for atom, energy in ra.items():
                    if atom not in ref_atoms:
                        ref_atoms[atom] = {}
                    ref_atoms[atom][bias] = energy * hartree_to_ev
        for atom, atom_data in ref_atoms.items():
            string += ',' + atom + ',,,' + ','.join([
                    '' if b not in atom_data else '%.3f'%atom_data[b] for b in biases]) + '\n'
        
        # add ref. molecules
        string += '\nRef. Mols.,,,,' + ','.join(biases) + '\n'
        for mol, bias_data in self.mols.items():
            energies = ['' if (b not in bias_data or not bias_data[b]['converged']) 
                        else '%.3f'%(bias_data[b]['final_energy']*hartree_to_ev) 
                        for b in biases]
            string += ','+mol+',,,'+','.join(energies) + '\n'
        
        # add surfaces 
        for surf, surf_data in analysis.items():
            mol_data = surf_data['analyzed']
            biases = self.order_bias(list(set([ b for mol in mol_data for b in mol_data[mol] ])))
            string += '\n' + surf + ',,,,' + ','.join(biases) + '\n'
            
            
            for mol, bias_data in mol_data.items():
                ii = 0
                mol_biases = list(set([b for b in bias_data]))
                try:
                    all_nsites = list(set([ n for b in mol_biases for n in bias_data[b]['all_site_data'] ]))
                except:
                    continue
                for nsite in all_nsites:
                    ads_site = [bias_data[b]['all_site_data'][nsite]['site'] for b in bias_data
                                if nsite in bias_data[b]['all_site_data']]
                    if not all([a == ads_site[0] for a in ads_site]):
                        print('Ads site changes with bias for '+mol+' at '+surf+' site '+nsite)
                        ads_site = 'var'
                    else:
                        ads_site = ads_site[0]
                        
                    bias_vals = ','.join([ '' if (b not in bias_data or nsite not in bias_data[b]['all_site_data'])
                                          else '%.3f'%(bias_data[b]['all_site_data'][nsite]['binding_energy']) 
                                          for b in biases])
                    string += ',' + (mol+',' if ii==0 else ',') + nsite + ',' + ads_site + ',' + bias_vals + '\n'
                    ii+=1
        
        # save data
        file = 'results/binding_energies.csv'
        with open(file,'w') as f:
            f.write(string)
    
    def system_analysis(self, data):
        '''
        Main function for analyzing converged data from scan_calcs function
        Functions:
            1) Creates analyzed.json file containing:
                - DONE: binding energies mapped over biases for each system (in eV)
                - DONE: get ads site (or sites, also get nsites)
                - DONE: get ads dist to site (in Ang)
                - DONE: get nelec delta (diff = ads_calc - surf - neutral_atoms)
                
                - Maybe: NEB barriers mapped over biases for each NEB system
        
        sys_analysis format:
            {ads_name: {bias: data}}
        '''
        sys_analysis = {}
        ads_data = data['adsorbed']
        surf_data = data['surf']
        
        for mol, bias_data in ads_data.items():
            mol_data = {}
            for bias, site_data in bias_data.items():
                min_site = None
                min_site_data = None
                min_energy = None
                all_site_data = {}
                for nsite, v in site_data.items():
                    surf_calc = surf_data[bias]
                    if not v['converged'] or not surf_calc['converged']:
                        continue # surf or ads not converged
                    binding_energy = self.system_dHf(v['final_energy'] - surf_calc['final_energy'], 
                                                     bias, self.formula_to_dic(mol))
                    if binding_energy is None:
                        continue
                    binding_site_data = self.get_ads_sites(Structure.from_dict(v['contcar']), mol, 
                                                           Structure.from_dict(surf_calc['contcar']), nsite)
                    if binding_site_data is None:
                        continue
                    binding_site = binding_site_data['ads_site']
                    nelec_diff = v['nfinal'] - surf_calc['nfinal'] - np.sum(
                                 [self.ref_atom_electrons[k]*v for k,v in self.formula_to_dic(mol).items()])
                    
                    all_site_data[nsite] = {'binding_energy': binding_energy, 'site': binding_site,
                                            #'site_data': binding_site_data, 
                                            'nelec_diff': nelec_diff, 
                                            'dist': binding_site_data['dist'],
                                            'nbonds': binding_site_data['nbonds'],
                                            'ads_atom': binding_site_data['ads']}
                    if min_energy is None or binding_energy < min_energy:
                        min_energy = binding_energy
                        min_site = (binding_site, nsite)
                        min_site_data = binding_site_data
                        min_nelec_diff = nelec_diff
                
                if min_energy is None:
                    mol_data[bias] = 'None'
                    continue
                mol_data[bias] = {'min_binding_energy': min_energy, 'energy_units': 'eV',
                                  'min_site_atom': min_site[0], 'min_site_id': min_site[1],
                                  'min_site_dist': min_site_data['dist'], 
                                  'min_site_nbonds': min_site_data['nbonds'],
                                  'ads_atom': min_site_data['ads'], 'nelec_diff': min_nelec_diff,
                                  'all_site_data': all_site_data}
                
            sys_analysis[mol] = mol_data
        return sys_analysis
    
    def get_ads_sites(self, ads_st, ads, surf_st, nsite, same_threshold = 1.5, max_bond = 3.3):
        # wrote this function without comments originally, no idea what's going on but it seems to work
        ads_dic = self.formula_to_dic(ads)
        site_dic = {}
        ads_list = [(k, i) for k,v in ads_dic.items() for i in range(v)]
        
        ads_sites = [site for site in ads_st.sites 
                     if not any([site.distance(s2) < same_threshold 
                     and site.species_string == s2.species_string for s2 in surf_st.sites])]
        
        for atom in ads_list:
            a = atom[0]
            i = atom[1]
            key = a+'_'+str(i)
            try:
                site = [s for s in ads_sites if s.species_string == a][i]
            except:
                print('Structure of adsorbate calc may be incorrect. Skipping.')
                return None
            site_dic[key] = {'site': site, 'dists': {str(ii).zfill(2): {'dist': site.distance(s2), 'site': s2}
                             for ii, s2 in enumerate(ads_st.sites) if s2 != site}, 'atom': a}
            site_dic[key]['bonds'] = {v['site'].species_string: v for k,v in site_dic[key]['dists'].items() 
                                      if v['dist'] < max_bond and v['site'] not in ads_sites}
            site_dic[key]['nbonds'] = len(site_dic[key]['bonds'])
            site_dic[key]['min_bond'] = {k: v for k,v in site_dic[key]['bonds'].items() if v['dist'] == min([
                                         b['dist'] for b in site_dic[key]['bonds'].values()])}
        
        min_bond = max_bond
        min_atom = ''
        for k1,v1 in site_dic.items():
            for k2,v2 in v1['min_bond'].items():
                if v2['dist'] < min_bond:
                    min_bond = v2['dist']
                    min_atom = v1['atom']
                    min_site = k2
                    min_sites = {k3: v3['dist'] for k3,v3 in v1['bonds'].items()}
                    min_nbonds = v1['nbonds']
        
        if min_atom == '':
            print('WARNING: No binding sites found: '+ads+' ('+nsite+')')
            return {'full_bond_dic': site_dic, 'ads': min_atom, 'ads_site': '', 'dist': 10, 'nbonds': 0}
        
        return {'full_bond_dic': site_dic, 'ads': min_atom, 'ads_site': min_site, 'all_sites': min_sites,
               'dist': min_bond, 'nbonds': min_nbonds}
    
    def ref_atoms(self):
        file = 'results/ref_atoms.json'
        if os.path.exists(file):
            with open(file,'r') as f:
                all_ref_atoms = json.load(f)
        else:
            all_ref_atoms = {'H': (['H3O','H2O'],[1,-1]), 
                         'N': (['N2'], [1/2]), 
                         'O': (['O2'], [1/2]),
                         'C': (['CO2','O2'], [1,-1]),}
        return all_ref_atoms
        
        
    def ref_energies(self, mol_dic, bias):
#        all_ref_atoms = {'H': (['H3O','H2O'],[1,-1]), 
#                         'N': (['N2'], [1/2]), 
#                         'O': (['O2'], [1/2]),
#                         'C': (['CO2','O2'], [1,-1]),}
        all_ref_atoms = self.ref_atoms()
        
        refs = {}
        for atom in mol_dic:
            if atom not in all_ref_atoms:
                print('Analysis Debug: atom references not known: '+atom)
                return None
            
            ref_mols = all_ref_atoms[atom][0]
            for mol in ref_mols:
                if mol not in self.mols:
                    print('Analysis Error: reference molecule not calculated: '+mol)
                    return None
            
            # get all biases shared by ref_mols
#            all_biases = [b for b in self.mols[ref_mols[0]] if 
#                          all([b in self.mols[r] and self.mols[r][b]['converged'] for r in ref_mols])]
            
            if not all([bias in self.mols[r] and self.mols[r][bias]['converged'] for r in ref_mols]):
                print('Analysis Error: not all needed biases run for mols '+str(ref_mols))
                return None
                        
            refs[atom] = (np.sum([all_ref_atoms[atom][1][i] 
                          * self.mols[mol][bias]['final_energy'] for i,mol in enumerate(ref_mols)]))
        return refs # in Hartree
    
    def system_dHf(self, total_energy, bias, mol_dic):
        refs = self.ref_energies(mol_dic, bias)
        if refs is None:
            return None
        return (total_energy - np.sum([refs[atom] * v for atom,v in mol_dic.items()])) * hartree_to_ev
    
    @property
    def ref_atom_electrons(self):
        return {'H': 1, 'O': 6, 'N': 5, 'C': 4, 'S': 6}
    
#    @property
    def reference_molecules(self, verbose = False):
        refs = {'H': {'refs': ['H2'], 'coeffs': [0.5]},
                'H2': {'refs': ['H2'], 'coeffs': [1]},
                'H2O': {'refs': ['H2O'], 'coeffs': [1]},
                'H3O':{'refs': ['H2O', 'H'], 'coeffs': [1, 1]}, 
                'H_H2O':{'refs': ['H2O', 'H'], 'coeffs': [1, 1]},
                'H_H3O':{'refs': ['H2O', 'H2'], 'coeffs': [1, 1]},
                'H2_H2O':{'refs': ['H2O', 'H2'], 'coeffs': [1, 1]},
                
                'O': {'refs': ['H2O','H2'], 'coeffs': [1,-1]},
                'OH': {'refs': ['O','H'], 'coeffs': [1,1]},
                'OOH': {'refs': ['O','H'], 'coeffs': [2,1]},
                'O2': {'refs': ['O2'], 'coeffs': [1]},
                
                'CO2': {'refs': ['CO2'], 'coeffs': [1]},
                'CO': {'refs': ['CO2','O'], 'coeffs': [1,-1]},
                'C': {'refs': ['CO2','O'], 'coeffs': [1,-2]},
                'CHO': {'refs': ['CO','H'], 'coeffs': [1,1]},
                'COH': {'refs': ['CHO'], 'coeffs': [1]},
                'COOH': {'refs': ['CO2','H'], 'coeffs': [1,1]},
                
                'OC': {'refs': ['CO'], 'coeffs': [1,]},
                'OCO': {'refs': ['CO2'], 'coeffs': [1]},
                'OCHO': {'refs': ['CO2','H'], 'coeffs': [1,1]},
                'HCOO': {'refs': ['CO2','H'], 'coeffs': [1,1]},
                'OCHOH': {'refs': ['CO2','H'], 'coeffs': [1,2]},
                'OCH': {'refs': ['CHO'], 'coeffs': [1]},
                'OCH2': {'refs': ['CO','H'], 'coeffs': [1,2]},
                'OCH3': {'refs': ['CO','H'], 'coeffs': [1,3]},
                'CHOH': {'refs': ['CO','H'], 'coeffs': [1,2]},
                'CH3OH': {'refs': ['CO','H'], 'coeffs': [1,2]},
                
                'N': {'refs': ['N2'], 'coeffs': [0.5]},
                'N2': {'refs': ['N2'], 'coeffs': [1]},
                'NN': {'refs': ['N2'], 'coeffs': [1]},
                'N2H': {'refs': ['N2','H'], 'coeffs': [1,1]},
                'NNH2': {'refs': ['N2','H'], 'coeffs': [1,2]},
                'NHNH': {'refs': ['N2','H'], 'coeffs': [1,2]},
                'NNH3': {'refs': ['N2','H'], 'coeffs': [1,3]},
                'NHNH2': {'refs': ['N2','H'], 'coeffs': [1,3]},
                'NH2NH2': {'refs': ['N2','H'], 'coeffs': [1,4]},
                'NH': {'refs': ['N2','H'], 'coeffs': [0.5,1]},
                'NH2': {'refs': ['N2','H'], 'coeffs': [0.5,2]},
                'NH3': {'refs': ['N2','H'], 'coeffs': [0.5,3]},
                
                'S2': {'refs': ['S8'], 'coeffs': [0.25]},
                'S4': {'refs': ['S8'], 'coeffs': [0.5]},
                'S6': {'refs': ['S8'], 'coeffs': [0.75]},
                'S8': {'refs': ['S8'], 'coeffs': [1]},
                }
        
        try:
            if os.path.exists('./molecules/refs.json'):
                with open('./molecules/refs.json', 'r') as f:
                    refs2 = json.load(f)
                
                for k,v in refs2.items():
                    refs[k] = v
        except:
            if verbose: print('Added refs not properly read from molecules/refs.json.')
        
        return refs
    
    def read_convergence(self, root):
        '''
        convergence example:
        step 1
        kpoints 1 1 1
        
        step 2
        kpoints 3 3 3
        '''
        with open(opj(root, 'convergence'),'r') as f:
            conv_txt = f.read()
        step = '1'
        add_step = False
        conv_dic = {}
        for line in conv_txt.split('\n'):
            if line == '' or line == ' ': 
                continue
#            print(line)
            if any(x in line for x in ['step ','Step ']):
                step = line.split()[1]
                if step == '0':
                    add_step = True
                if add_step:
                    step  = str(int(step)+1) # update so steps are always indexed to 1
                if step not in conv_dic:
                    conv_dic[step] = {}
            else:
                cmd, val = line.split()[0], ' '.join(line.split()[1:])
                if cmd in conv_dic[step]:
                    # multiple versions of same command (i.e. pdos)
                    if type(conv_dic[step][cmd]) == list:
                        conv_dic[step][cmd].append(val)
                    else:
                        conv_dic[step][cmd] = [conv_dic[step][cmd], val]
                else:
                    # single version of command
                    conv_dic[step][cmd] = val
        return conv_dic
    
    def write_convergence(self, root, conv_dic):
        txt = ''
        assert '1' in conv_dic, 'METAERROR: Value Step 1 not in convergence dictionary, cannot write file.'
        for step, step_dic in conv_dic.items():
            if step == '1':
                txt += 'step 1\n'
            else:
                txt += '\nstep '+step + '\n'
            
            for cmd, val in step_dic.items():
                if type(val) == list:
                    for v in val:
                        txt += cmd + ' ' + v + '\n'
                else:
                    txt += cmd + ' ' + val + '\n'
        
        # write convergence file
        with open(opj(root, 'convergence'), 'w') as f:
            f.write(txt)
            
    def read_eigStats(self, folder):
        file = opj(folder, 'eigStats')
        if ope(file):
            with open(file, 'r') as f:
                eig = f.read()
            eigs = {'units': 'H'}
            for line in eig.split('\n'):
                if line in ['',' ']:
                    continue
                key = line.split(': ')[0]
                val = float(line.split(': ')[1].split(' ')[0])
                eigs[key] = val
            return eigs
        return {}
    
    def read_forces(self, folder):
        file = opj(folder, 'forces')
        if ope(file):
            with open(file, 'r') as f:
                forces = f.read()
            f = []
            for ii,line in enumerate(forces.split('\n')):
                if ii == 0 or len(line) == 0: continue
                atom = line.split()[1]
                sd = int(line.split()[5])
                x = float(line.split()[2])
                y = float(line.split()[3])
                z = float(line.split()[4])
                f.append({'atom': atom, 'number': ii, 'SD': sd, 'x': x, 'y': y, 'z': z})
            return f
        return []
    
    def get_jdos(self, folder, plot = False):
        cwd = os.getcwd()
        os.chdir(folder)
        
        if os.path.exists('jpdos.json'):
            with open('jpdos.json','r') as f:
                jdos = json.load(f)
            os.chdir(cwd)
            return jdos
        try:
            subprocess.call('get_jdos.py', shell=True)
            with open('jpdos.json','r') as f:
                jdos = json.load(f)
        except:
            print('ERROR: Cannot read DOS info for:', folder)
            jdos = {}
        if plot:
            print('METAERROR: jdos plotting not yet added! Contact Nick.')
            #self.plot_jdos(jdos)
        os.chdir(cwd)
        return jdos
        
    def formula_to_dic(self, formula):
        assert '.' not in formula, 'ERROR: Cannot include decimals in formula. Only integers allowed.'
        if '(' in formula:
            first = formula.split('(')[0]
            second = ''.join(formula.split('(')[1:]).split(')')[0]
            coeff = ''.join(formula.split('(')[1:]).split(')')[1]
            d = self.formula_to_dic(first)
            for k, v in self.formula_to_dic(second).items():
                if k in d:
                    d[k] += v*int(coeff)
                else:
                    d[k] = v*int(coeff)
            return d
        
        # split formula to list by number
        spl = [s for s in re.split('(\d+)',formula) if s != '']
        els_list = []; el_frac = [];
        for el in spl: 
            if sum(1 for c in el if c.isupper()) == 1 and len(el) < 3:
                # I am a single element, add me to the list please
                els_list.append(el)
                el_frac.append(1)
            elif el.isdigit() == True:
                # I am a number and should be calculated in to self.el_frac
                el_frac[-1] = el_frac[-1] * int(el)
            elif sum(1 for c in el if c.isupper()) > 1:
                # I am more than one element and am here to cause trouble
                els_split = re.sub( r"([A-Z])", r" \1", el).split()
                els_list += els_split
                for i in els_split:
                    el_frac.append(1)
        # combine duplicate elements
        for el in els_list:
            inds = [i for i, j in enumerate(els_list) if j == el]
            if len(inds) > 1:
                drop_ind = inds[1:]
                drop_ind = drop_ind[::-1]
                el_frac[inds[0]] = sum(el_frac[inds[i]] for i in range(len(inds)))
                for d in drop_ind: del els_list[d]; del el_frac[d]
        assert len(els_list) == len(el_frac), 'ERROR: elements and coefficients not balanced in "formula_to_dic".'
#        assert all([ k in self.el_masses.keys() for k in els_list ]), 'ERROR: unknown element detected: '+formula
        return {els_list[i]: el_frac[i] for i,_ in enumerate(els_list)}
        

    
    
