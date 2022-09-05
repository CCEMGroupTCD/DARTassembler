
class Employee:

    raise_amt = 1.04

    def __init__(self, first, last, pay):
        self.first = first
        self.last = last
        self.email = first + '.' + last + '@email.com'
        self.pay = pay

    def fullname(self):
        return '{} {}'.format(self.first, self.last)

    def apply_raise(self):
        self.pay = int(self.pay * self.raise_amt)


class Employee2(Employee):
    def print_(self):
        print(self.raise_amt)


if __name__ == '__main__':
    dev_2 = Employee('Test', 'Employee', 60000)
    dev_1 = Employee2('Test', 'Employee', 60000)
    print("done")