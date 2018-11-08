class MathOperations:
    @staticmethod
    def testAddition (x, y):
        return x + y
    @staticmethod
    def testMultiplication (a, b):
        return a * b
    @staticmethod
    def multadd(n,m):
        A = MathOperations.testAddition(n,m)
        B = MathOperations.testAddition(n, -m)
        return MathOperations.testMultiplication(A,B)
    
        
tmp = MathOperations
print tmp.multadd(3,10)