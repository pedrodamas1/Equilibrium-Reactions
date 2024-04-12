
class Atom:
    
    def __init__(self, name: str) -> None:
        self.name = name

    def __str__(self) -> str:
        return self.name
    
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.name}')"


if __name__ == '__main__':
    iron = Atom('Fe')
    print(str(iron))
    print(repr(iron))
