class StateLiteral:
    def __init__(self, state: int):
        self.state = state

    def __repr__(self):
        return f'[{self.state}]'

    def __eq__(self, other):
        if isinstance(other, StateLiteral):
            return self.state == other.state
        return False

    def __hash__(self):
        return hash(self.state)

ZERO_STATE = StateLiteral(0)
NUNITY_STATE = StateLiteral(-1)