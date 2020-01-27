import numpy as np
import math
import pprint
import argparse

def lattice_vector_size(array):
	return math.sqrt(array[0]**2+array[1]**2+array[2]**2)
def vector_distance(array_1,array_2):
	return lattice_vector_size(array_1-array_2)
def cartesian_transformation(lattice_vector_array,atom_coordinate_array):
	temp_cartesian=[]
	for i in range(3):
		temp_cartesian.append(lattice_vector_array[0][i]*atom_coordinate_array[0]+
			lattice_vector_array[1][i]*atom_coordinate_array[1]+
			lattice_vector_array[2][i]*atom_coordinate_array[2])
	temp_cartesian=np.array(temp_cartesian,dtype=np.float64)
	return temp_cartesian
def read_poscar(open_file_line,line_number): #x는 오픈한 것, y는 라인 넘버
	temp=open_file_line[line_number-1]
	temp=temp.replace('\n','') #n을 날린다
	temp=temp.split() #공백을 기준으로 리스트 만든거임
	return temp
def read_poscar_num(open_file_line,line_number):
	temp=open_file_line[line_number-1]
	temp=temp.replace('\n','') #n을 날린다
	temp=temp.split()
	for i in range(0,3):
		if float(temp[i])<0:
			temp[i]=float(temp[i])+1
	temp=np.array(temp,dtype=np.float64)
	list(temp)
	return temp
def lattice_vector(open_file_line):
	temp=[]
	for i in range(2,5):
		line=open_file_line[i].replace('\n','')
		line=line.split()
		temp.append(line)
	temp=np.array(temp,dtype=np.float64)
	return temp
def atom_coor_direc(open_file_line,atom_list,atom_num_list):
	total_line_coor=[]
	for i in range(9,len(open_file_line)+1):
		total_line_coor.append(read_poscar_num(open_file_line,i))
	atom_list_with_num_line={}
	for i in range(0,len(atom_list)):
		new_coor=np.array(total_line_coor[0:atom_num_list[i]],dtype=np.float64)
		atom_list_with_num_line['%s'%atom_list[i]]=new_coor
		total_line_coor=total_line_coor[atom_num_list[i]:]
	return atom_list_with_num_line
def atom_coor_cartesian(open_file_line,atom_list,atom_num_list,lattice_vector_):
	temp=atom_coor_direc(open_file_line,atom_list,atom_num_list)
	atom_list_with_num_line={}
	for i in temp:
		temp_atom_list=[]
		for j in temp[i]:
			j=cartesian_transformation(lattice_vector_,j)
			temp_atom_list.append(j)
		temp_atom_list_1=[]
		for k in temp_atom_list:
			temp_atom_list_1.append(k)
		atom_list_with_num_line[i]=np.array(temp_atom_list_1,dtype=np.float64)
	return atom_list_with_num_line
def super_atom_coor_direc(open_file_line,atom_list,atom_num_list,shift_list):
	temp=atom_coor_direc(open_file_line,atom_list,atom_num_list)
	atom_list_with_num_line={}
	for i in temp:
		temp_atom_list=[]
		for j in temp[i]:
			for k in shift_list:
				temp_atom_list.append(j+k)
		atom_list_with_num_line[i]=np.array(temp_atom_list,dtype=np.float64)
	return atom_list_with_num_line
def super_atom_coor_cartesian(open_file_line,atom_list,atom_num_list,shift_list,lattice_vector_):
	temp=super_atom_coor_direc(open_file_line,atom_list,atom_num_list,shift_list)
	atom_list_with_num_line={}
	for i in temp:
		temp_atom_list=[]
		for j in temp[i]:
			j=cartesian_transformation(lattice_vector_,j)
			temp_atom_list.append(j)
		temp_atom_list_1=[]
		for k in temp_atom_list:
			temp_atom_list_1.append(k)
		atom_list_with_num_line[i]=np.array(temp_atom_list_1,dtype=np.float64)
	return atom_list_with_num_line
def find_poly_for_original(lines,atom_list,atom_num_list,latt_vec,center_atom,vertex_atom,max_bond_length):
	shift_list=([[x,y,z] for x in range(-1,2) for y in range(-1,2) for z in range(-1,2)])
	shift_list=np.array(shift_list)
	original_center=[]
	super_vertex=[]
	for atom in atom_coor_cartesian(lines,atom_list,atom_num_list,latt_vec)[center_atom]:
		original_center.append(atom)
	for atom in super_atom_coor_cartesian(lines,atom_list,atom_num_list,shift_list,latt_vec)[vertex_atom]:
		super_vertex.append(atom)
	center_vertex=[]
	for atom_c in original_center:
		vertex_atoms=[]
		for atom_v in super_vertex:
			if lattice_vector_size(atom_c-atom_v) < max_bond_length:
				vertex_atoms.append([atom_v,lattice_vector_size(atom_c-atom_v)])
		center_vertex.append([atom_c,vertex_atoms])
	return center_vertex
def find_poly_for_super(lines,atom_list,atom_num_list,latt_vec,center_atom,vertex_atom,max_bond_length):
	shift_list_c=([[x,y,z] for x in range(-1,2) for y in range(-1,2) for z in range(-1,2)])
	shift_list_c=np.array(shift_list_c)
	shift_list_v=([[x,y,z] for x in range(-2,3) for y in range(-2,3) for z in range(-2,3)])
	shift_list_v=np.array(shift_list_v)
	super_center=[]
	super_vertex=[]
	for atom in super_atom_coor_cartesian(lines,atom_list,atom_num_list,shift_list_c,latt_vec)[center_atom]:
		super_center.append(atom)
	for atom in super_atom_coor_cartesian(lines,atom_list,atom_num_list,shift_list_c,latt_vec)[vertex_atom]:
		super_vertex.append(atom)
	super_center_vertex=[]
	for atom_c in super_center:
		vertex_atoms=[]
		for atom_v in super_vertex:
			if lattice_vector_size(atom_c-atom_v) < max_bond_length:
				vertex_atoms.append([atom_v,lattice_vector_size(atom_c-atom_v)])
		super_center_vertex.append([atom_c,vertex_atoms])
	return super_center_vertex
def calculating_off_centering(center_vertex,center_atom):
	total_off_centering=0
	poly_num=0
	for octahedral in center_vertex:
		poly_num+=1
		weight_center=np.array([0,0,0],dtype=np.float64)
		mean_bond_length=0
		for vertex_atom_num in range(len(octahedral[1])):
			#print(octahedral[0],octahedral[1][vertex_atom_num][0],octahedral[1][vertex_atom_num][1])
			weight_center+=octahedral[1][vertex_atom_num][0]
			mean_bond_length+=octahedral[1][vertex_atom_num][1]
		weight_center=weight_center/len(octahedral[1])
		mean_bond_length=mean_bond_length/len(octahedral[1])
		print('Poly%s\n%s           : %s %s %s'%(poly_num,center_atom,octahedral[0][0],octahedral[0][1],octahedral[0][2]))
		print('weight_center: %s %s %s'%(weight_center[0],weight_center[1],weight_center[2]))
		off_centering=vector_distance(octahedral[0],weight_center)*100/mean_bond_length
		total_off_centering+=off_centering
		vector_minus=octahedral[0]-weight_center
		print('vector       : %s %s %s'%(vector_minus[0],vector_minus[1],vector_minus[2]))
		print('Off centering: %0.3f (%%)\n'%off_centering)
	total_off_centering=total_off_centering/len(center_vertex)
	print('Total off centering: %0.3f (%%)\n'%total_off_centering)

def center_atom_lattice(center_vertex,super_center_vertex):
	distance_between_center=[]
	for original_atom in center_vertex:
		distance_between_center_temp=[]
		for super_atom in super_center_vertex:
			if np.array_equal(original_atom[0],super_atom[0]):
				pass
			else:
				distance_between_center_temp.append([super_atom[0],lattice_vector_size(original_atom[0]-super_atom[0])])
		distance_between_center_temp=sorted(distance_between_center_temp,key=lambda x: x[1])
		distance_between_center.append([original_atom[0],distance_between_center_temp])
	pprint.pprint(distance_between_center)

def print_bondlength(center_vertex,center_atom,vertex_atom):
	poly_num=0
	for polyhedra in center_vertex:
		poly_num+=1
		print("Poly %i, %s - %s bondlengths\n"%(poly_num,center_atom,vertex_atom))
		for vertex_atoms in polyhedra[1]:
			print(vertex_atoms[1])
		print('')
def calculating_offcenter_vector_sum(center_vertex,center_atom):
	total_off_centering=0
	poly_num=0
	vector_sum=[]
	bondlengths_sum=[]
	for octahedral in center_vertex:
		poly_num+=1
		weight_center=np.array([0,0,0],dtype=np.float64)
		bond_length=0
		for vertex_atom_num in range(len(octahedral[1])):
			#print(octahedral[0],octahedral[1][vertex_atom_num][0],octahedral[1][vertex_atom_num][1])
			weight_center+=octahedral[1][vertex_atom_num][0]
			bond_length+=octahedral[1][vertex_atom_num][1]
		weight_center=weight_center/len(octahedral[1])
		bondlengths_sum.append(bond_length)
		vector_sum.append((octahedral[0]-weight_center))
	mean_bond_length=sum(bondlengths_sum)/(poly_num*len(octahedral[1]))
	vector_sum=sum(vector_sum)
	print('vector sum: %s %s %s'%(vector_sum[0],vector_sum[1],vector_sum[2]))
	print('mean bond length: ',mean_bond_length)
def main(args):
	center_atom=args.center_atom
	vertex_atom=args.vertex_atom
	max_bond_length=args.max_bond_length
	file_name=args.file_name
	with open("%s"%file_name,"r") as f:
		lines=f.readlines()
	atom_num_list=list(np.array(read_poscar(lines,7),dtype=np.int))
	lines=lines[:8+sum(atom_num_list)]
	atom_list=read_poscar(lines,6)
	latt_vec=lattice_vector(lines)
	latt_par=[]
	for i in range(3):
		latt_par.append(lattice_vector_size(latt_vec[i]))
	latt_par=np.array(latt_par,dtype=np.float64)
	latt_par=latt_par.reshape(3,1)
	direct_coor=atom_coor_direc(lines,atom_list,atom_num_list)
	cartesian_coor=atom_coor_cartesian(lines,atom_list,atom_num_list,latt_vec)
	center_vertex=find_poly_for_original(lines,atom_list,atom_num_list,latt_vec,center_atom,vertex_atom,max_bond_length)
	# print bondlengths
	#print_bondlength(center_vertex,center_atom,vertex_atom)
	# off centering calculation
	print('')
	calculating_off_centering(center_vertex,center_atom)
	calculating_offcenter_vector_sum(center_vertex,center_atom)
	#super_center_vertex=find_poly_for_super(lines,atom_list,atom_num_list,latt_vec,center_atom,vertex_atom,max_bond_length)
	#center_atom_lattice(center_vertex,super_center_vertex)

parser = argparse.ArgumentParser(description='devide tasks')
parser.add_argument("-c","-center",type=str,dest="center_atom",
	help='Polyhedra center atom name')
parser.add_argument("-v","-vertex",type=str,dest="vertex_atom",
	help='Polyhedra vertex atom name')
parser.add_argument("-l","-length",type=float,default=3,
	dest="max_bond_length",help="Set maximum bond length")
parser.add_argument("-f","-file_name",type=str,dest="file_name",
	default='POSCAR',help='Enter POSCAR file name')
parser.set_defaults(func=main)
args=parser.parse_args()
args.func(args)