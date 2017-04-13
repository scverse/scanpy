
def get_M_row(i, evals, rbasis, lbasis):
    m_i = 0 
    for l in range(1, evals.size):
        m_i += evals[l]/(1-evals[l]) * rbasis[i, l] * lbasis[:, l]
    m_i += rbasis[i, 0] * lbasis[:, 0]
    return m_i
