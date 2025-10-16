import math
import pathlib


def load_system(path: pathlib.Path) -> tuple[list[list[float]], list[float]]:
    system=path.read_text()
    lines=system.splitlines()
    free_terms=[]
    system_matrix=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    for index,line in enumerate(lines):
        split_equal=line.split('=')
        split_equal[0]=split_equal[0]
        free_terms.append(float(split_equal[1]))
        left_part=split_equal[0].replace(' ','')
        split_x=left_part.split('x')
        if len(split_x) == 1:
            system_matrix[index][0] = 0.0
            split_y=split_x[0].split('y')
        else:
            system_matrix[index][0]=1.0 if split_x[0] in ('','+') else (-1.0 if split_x[0]=='-' else float(split_x[0]))
            split_y=split_x[1].split('y')
        if len(split_y) == 1:
            system_matrix[index][1] = 0.0
            split_z=split_y[0].split('z')
        else:
            system_matrix[index][1]=1.0 if split_y[0] in ('','+') else (-1.0 if split_y[0]=='-' else float(split_y[0]))
            split_z=split_y[1].split('z')
        if len(split_z) == 1:
            system_matrix[index][2] = 0.0
        else:
            system_matrix[index][2] = 1.0 if split_z[0] in ('', '+') else (-1.0 if split_z[0] == '-' else float(split_z[0]))
    return system_matrix,free_terms

def determinant(matrix: list[list[float]]) -> float:
    a11=matrix[0][0]*(matrix[1][1]*matrix[2][2]-matrix[1][2]*matrix[2][1])
    a12=matrix[0][1]*(matrix[1][0]*matrix[2][2]-matrix[1][2]*matrix[2][0])
    a13=matrix[0][2]*(matrix[1][0]*matrix[2][1]-matrix[1][1]*matrix[2][0])
    return a11-a12+a13

def trace(matrix: list[list[float]]) -> float:
    return sum(matrix[i][i] for i in range(len(matrix)))


def norm(vector: list[float]) -> float:
    return math.sqrt(sum(x * x for x in vector))


def transpose(matrix: list[list[float]]) -> list[list[float]]:
    return [[ matrix[row][col] for row in range(len(matrix)) ] for col in range(len(matrix[0]))]

def multiply(matrix: list[list[float]], vector: list[float]) -> list[float]:
    product=[0.0 for _ in range(len(vector))]
    for row in range(len(matrix)):
        for col in range(len(matrix[row])):
            product[row]+=matrix[row][col]*vector[col]
    return product

def solve_cramer(matrix: list[list[float]],vector: list[float]) -> list[float]:
   solution=[]
   delta=determinant(matrix)
   if delta==0.0:
       raise ValueError("The system is not uniquely solvable system")
   for i in range(3):
    for index in range(len(vector)):
        matrix[index][i],vector[index]=vector[index],matrix[index][i]
    solution.append(determinant(matrix)/delta)
    for index in range(len(vector)):
        matrix[index][i], vector[index] = vector[index], matrix[index][i]

   return solution

def minor(matrix: list[list[float]], i: int, j: int) -> list[list[float]]:
   minor_matrix=[]
   for row_idx,row in enumerate(matrix):
       if row_idx==i:
           continue
       minor_row=[]
       for col_idx,col in enumerate(row):
           if col_idx==j:
               continue
           minor_row.append(col)
       minor_matrix.append(minor_row)
   return minor_matrix

def cofactor(matrix: list[list[float]]) -> list[list[float]]:
    n=len(matrix)
    cofactor_matrix=[]
    for i in range(n):
        cofactor_row=[]
        for j in range(n):
            minor_matrix=minor(matrix,i,j)
            delta_minor=(-1)**(i+j)
            delta_minor*=(minor_matrix[0][0]*minor_matrix[1][1]-minor_matrix[0][1]*minor_matrix[1][0])
            cofactor_row.append(delta_minor)
        cofactor_matrix.append(cofactor_row)
    return cofactor_matrix


def adjoint(matrix: list[list[float]]) -> list[list[float]]:
     return cofactor(transpose(matrix))

def solve(matrix: list[list[float]], vector: list[float]) -> list[float]:
    solution=multiply(adjoint(matrix),vector)
    delta=determinant(matrix)
    return [val/delta for val in solution]

A,B=load_system(pathlib.Path("/Labs/system.txt"))

print(f"{A=} {B=}")

print(f"{determinant(A)=}")

print(f"{trace(A)=}")

print(f"{norm(B)=}")

print(f"{transpose(A)=}")

print(f"{multiply(A, B)=}")

print(f"{solve_cramer(A, B)=}")

print(f"{solve(A, B)=}")




