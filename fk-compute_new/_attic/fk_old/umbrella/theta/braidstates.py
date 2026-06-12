import symengine
import sympy

class BraidStates:
    def __init__(self, braid, inversion_data=None):
        self.braid = braid
        self.max_strand = max(abs(x) for x in braid)
        self.strands = list(range(0,self.max_strand+1))
        self.n_strands = self.max_strand + 1
        self.n_crossings = len(braid)
        self.braid_group_generators = list(range(1,self.max_strand+1))
        self.crossing_signs = []
        for g in braid:
            if g < 0:
                self.crossing_signs.append(-1)
            else:
                self.crossing_signs.append(+1)
        self.writhe = sum(self.crossing_signs)

        self.top_input_state_locations = [(abs(braid[index]) - 1, index) for index in range(self.n_crossings)]
        self.bottom_input_state_locations = [(x + 1, y) for (x, y) in self.top_input_state_locations]
        self.loc = [[[(x, y + 1), (x + 1, y + 1)], [(x, y), (x + 1, y)]] for (x, y) in self.top_input_state_locations]

        self.T = symengine.Symbol('T')
        self.T1 = symengine.Symbol('T1')
        self.T2 = symengine.Symbol('T2')

        current_list = []
        self.strand_locations = []
        current_loc = [0, 0]
        done = False
        while not done:
            current_list.append(tuple(current_loc))
            if tuple(current_loc) in self.bottom_input_state_locations:
                self.strand_locations.append(current_list)
                current_list = []
                current_loc = [current_loc[0] - 1, current_loc[1] + 1]
            elif tuple(current_loc) in self.top_input_state_locations:
                self.strand_locations.append(current_list)
                current_list = []
                current_loc = [current_loc[0] + 1, current_loc[1] + 1]
            elif current_loc[1] == self.n_crossings:
                current_loc = [current_loc[0], 0]
            else:
                current_loc = [current_loc[0], current_loc[1] + 1]
            if current_loc == [0, 0]:
                self.strand_locations.append(current_list)
                done = True
        self.position_strands = dict()
        for index in range(len(self.strand_locations)):
            for index_ in range(len(self.strand_locations[index])):
                self.position_strands[self.strand_locations[index][index_]] = index
        self.strand_endpoints = [[x[0], x[-1]] for x in self.strand_locations]
        self.endpoint_crossing_indices = [[x[0][1] - 1, x[-1][1]] for x in self.strand_endpoints]
        self.n_s = len(self.strand_locations)

    def compute_A(self):
        self.A = symengine.eye(self.n_s)
        for index in range(self.n_crossings):
            s = self.crossing_signs[index]
            if s == 1:
                i = self.position_strands[self.loc[index][1][0]]
                j = self.position_strands[self.loc[index][1][1]]
                i_ = self.position_strands[self.loc[index][0][1]]
                j_ = self.position_strands[self.loc[index][0][0]]
            else:
                i = self.position_strands[self.loc[index][1][1]]
                j = self.position_strands[self.loc[index][1][0]]
                i_ = self.position_strands[self.loc[index][0][0]]
                j_ = self.position_strands[self.loc[index][0][1]]
            self.A[i * self.n_s + i_] -= self.T ** (s)
            self.A[i * self.n_s + j_] += self.T ** (s) - 1
            self.A[j * self.n_s + j_] -= 1

    def compute_G(self):
        self.G = self.A.inv()
        for index in range(self.n_s * self.n_s):
            self.G[index] = self.G[index]

    def compute_Alexander(self):
        T = symengine.Symbol('T')
        self.Alexander = T**((self.n_strands - 1 - self.writhe) / 2) * self.A.det()

    def F1(self):
        summand = 0
        for index in range(self.n_crossings):
            if self.crossing_signs[index] == 1:
                i = self.position_strands[self.loc[index][1][0]]
                j = self.position_strands[self.loc[index][1][1]]
                s = 1
            else:
                j = self.position_strands[self.loc[index][1][0]]
                i = self.position_strands[self.loc[index][1][1]]
                s = -1
            gii = [0, self.G[i * self.n_s + i].subs({self.T:self.T1}), self.G[i * self.n_s + i].subs({self.T:self.T2}), self.G[i * self.n_s + i].subs({self.T:self.T1*self.T2})]
            gij = [0, self.G[i * self.n_s + j].subs({self.T:self.T1}), self.G[i * self.n_s + j].subs({self.T:self.T2}), self.G[i * self.n_s + j].subs({self.T:self.T1*self.T2})]
            gji = [0, self.G[j * self.n_s + i].subs({self.T:self.T1}), self.G[j * self.n_s + i].subs({self.T:self.T2}), self.G[j * self.n_s + i].subs({self.T:self.T1*self.T2})]
            gjj = [0, self.G[j * self.n_s + j].subs({self.T:self.T1}), self.G[j * self.n_s + j].subs({self.T:self.T2}), self.G[j * self.n_s + j].subs({self.T:self.T1*self.T2})]
            T1 = self.T1 ** s
            T2 = self.T2 ** s
            summand += s*(1/2 - gii[3] + T2*gii[1]*gji[2] - T2*gjj[3]*gji[2] - (T2 - 1)*gii[3]*gji[2] + (T1*T2 - 1)*gji[2]*gji[3] - gii[1]*gjj[2] + 2*gii[3]*gjj[2] + gii[1]*gjj[3] - gii[2]*gjj[3]) + s/(T2 - 1)*((T1 - 1)*T2*(gjj[3]*gji[1] - gjj[2]*gji[1] + T2*gji[1]*gji[2]) + (T1*T2 - 1)*(gji[3] - T2*gii[1]*gji[3] + gij[2]*gji[3] + (T2 - 2)*gjj[2]*gji[3]) - (T1 - 1)*(T2 + 1)*(T1*T2 - 1)*gji[1]*gji[3])
        return summand 
    
    def F2(self):
        summand = 0
        for index0 in range(self.n_crossings):
            if self.crossing_signs[index0] == 1:
                i0 = self.position_strands[self.loc[index0][1][0]]
                j0 = self.position_strands[self.loc[index0][1][1]]
                s0 = 1
            else:
                j0 = self.position_strands[self.loc[index0][1][0]]
                i0 = self.position_strands[self.loc[index0][1][1]]
                s0 = -1
            for index1 in range(self.n_crossings):
                if self.crossing_signs[index1] == 1:
                    i1 = self.position_strands[self.loc[index1][1][0]]
                    j1 = self.position_strands[self.loc[index1][1][1]]
                    s1 = 1
                else:
                    j1 = self.position_strands[self.loc[index1][1][0]]
                    i1 = self.position_strands[self.loc[index1][1][1]]
                    s1 = -1
                g1j1i0 = self.G[j1 * self.n_s + i0].subs({self.T:self.T1})
                gj0i1 = self.G[j0 * self.n_s + i1].subs({self.T:self.T1*self.T2})
                gi1i0 = self.G[i1 * self.n_s + i0].subs({self.T:self.T2})
                gj1j0 = self.G[j1 * self.n_s + j0].subs({self.T:self.T2})
                g2j1i0 = self.G[j1 * self.n_s + i0].subs({self.T:self.T2})
                gi1j0 = self.G[i1 * self.n_s + j0].subs({self.T:self.T2})
                summand += s1*(self.T1 ** (s0) - 1)*((self.T1*self.T2) ** (s1) - 1)*g1j1i0*gj0i1/(self.T2 ** (s1) - 1)*(self.T2 ** (s0) * gi1i0 + gj1j0 - self.T2 ** (s0) * g2j1i0 - gi1j0)
        return summand
                
    def F3(self):
        summand = 0
        for index in range(self.n_s):
            if self.strand_endpoints[index][0][1] > self.strand_endpoints[index][1][1] and self.strand_endpoints[index][0][0] != 0:
                summand += -1 * (self.G[index * self.n_s + index].subs({self.T:self.T1*self.T2}) - 1/2)
        return summand
    
    def compute_Theta(self):
        self.Theta = self.Alexander.subs({self.T:self.T1}) * self.Alexander.subs({self.T:self.T2}) * self.Alexander.subs({self.T:self.T1*self.T2}) * (self.F1() + self.F2() + self.F3())

if __name__ == "__main__":

    # braid_states = BraidStates([1,1,-2,-1,-1,3,2,-1,2,2,3])
    braid_states = BraidStates([1,1,1])
    braid_states.compute_A()
    braid_states.compute_G()
    braid_states.compute_Alexander()
    print(braid_states.Alexander.simplify())
    braid_states.compute_Theta()
    print(braid_states.Theta.simplify())
    print('\n')
    print(sympy.cancel(braid_states.Theta.simplify().subs({braid_states.T1:1})))