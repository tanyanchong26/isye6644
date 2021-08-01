class random:
    
    def _init_(self):
        self.random_number = []
        self.bins = []
    
    #%% Generators
    
    def generator_desert_island(self, seed, quantity):    
        """
        The first paramater is the initial seed. It must be an integer.
        The second parameter is the size of the vector of random numbers desired.
        """    
        if seed < 1:
            raise AssertionError("Seed must be larger than or equal to 1!")
    
        if not isinstance(seed,int):
            raise AssertionError("Seed must be an integer!")
        
        if not isinstance(quantity,int):
            raise AssertionError("Vector Length must be an integer!")
    
        if quantity < 1:
            raise AssertionError("Vector Length must be greater than or equal to 1!")
    
        X = []
    
        for i in range(0, quantity+1):
            if len(X)==0:
                X.append((16807*seed)%(2**31-1))
            else:
                X.append((16807*X[-1])%(2**31-1))    
       
        self.random_number = [i/(2**31-1) for i in X[1:]]
    
    def generator_randu(self, seed, quantity):
    
        """
        The first paramater is the initial seed. It must be an integer.
        The second parameter is the size of the vector of random numbers desired.
        """
    
        if seed < 1:
            raise AssertionError("Seed must be larger than or equal to 1!")
    
        if not isinstance(seed,int):
            raise AssertionError("Seed must be an integer!")
        
        if not isinstance(quantity,int):
            raise AssertionError("Vector Length must be an integer!")
    
        if quantity < 1:
            raise AssertionError("Vector Length must be greater than or equal to 1!")
    
        X = []
    
        for i in range(0, quantity):
            if len(X)==0:
                X.append((65539*seed)%(2**31))
            else:
                X.append((65539*X[-1])%(2**31))    
       
        self.random_number = [i/(2**31) for i in X]
    
    def generator_tausworthe(self, quantity):    
        """
        Only accepts 1 parameter, which is the number of random numbers desired.
        """
    
        if not isinstance(quantity, int) or quantity<=0:
            raise AssertionError("Quantity must be greater than 0 and must be an integer value!")
    
        U = []
        r = 2
        q = 35
        l = 35
        n = quantity
    
        B = [True]*q

        len_required = l*n

        while len(B)<len_required:
            B.append(B[-r] != B[-q])

        for i in range(0,len(B),l):
            result = B[i:(i+l)]
            result_inv = result[::-1]
            seq_number = [i for i in range(0, l)]
            random_number = sum([i*(2**j) for i,j in zip(result_inv, seq_number)])/2**l
            U.append(random_number)
    
        self.random_number = U
    
    def generator_lecuyer(self, seed_1, seed_2, quantity):    
        """
        The two parents are initial seeds. They must be in list with length of 3 and 
        each element must be larger than 0 in each of the list of seeds.
        
        The size parameter indicates the length of random numbers required. It must be an integer value
        greater than 1.
        """
    
        if not isinstance(seed_1,list):
            raise AssertionError("Seed 1 has to be a List!")
    
        if not isinstance(seed_2,list):
            raise AssertionError("Seed 2 has to be a List!")
    
        if len(seed_1) != 3:
            raise AssertionError("Length of Seed 1 must be 3!")
    
        if len(seed_2) != 3:
            raise AssertionError("Length of Seed 2 must be 3!")
    
        for i in seed_1:
            if i <= 0:
               raise AssertionError("Seeds less than or equal to 0 found in Seed 1!")
    
        for i in seed_2:
            if i <= 0:
                raise AssertionError("Seeds less than or equal to 0 found in Seed 2!")
    
        if not isinstance(quantity, int) or quantity<=0:
            raise AssertionError("Quantity must be greater than 0 and must be an integer value!")
    
        X1 = seed_1.copy()
        X2 = seed_2.copy()
    
        n = 0
        random_numbers = []
    
        while n<(quantity+1):        
            X1i = (1403580*X1[-2]-810728*X1[-3])%(2**32-209)
            X1.append(X1i)
    
            X2i = (527612*X2[-1]-1370589*X2[-3])%(2**32-22853)
            X2.append(X2i)
    
            Yi = (X1i - X2i)%(2**32-209)
            Ri = Yi/(2**32-209)
            random_numbers.append(Ri)
            n += 1
    
        self.random_number = random_numbers[1:]
    
    # Auxiliary Functions
    
    def generator_number_binning(self, segments = 30):
        """
        The first parameter is a list of randomly generated numbers from 0 to 1.
        The second parameter is the number of bins desired to be segmented into.
        """ 
        
        import numpy as np
        
        if not isinstance(self.random_number,list):
            raise AssertionError("The input must be a list!")
    
        if segments <= 1 or not isinstance(segments, int):
            raise AssertionError("The number of bins must be an integer!")
    
        R = np.array(self.random_number)
        bin_results = {}
        bins = np.linspace(0,1,segments+1)
        segment_list = np.digitize(R, bins)
    
        for i in segment_list:
            if i not in bin_results:
                bin_results[i] = 1
            else:
                bin_results[i] += 1
    
        bin_results = dict(sorted(bin_results.items()))
    
        self.bins = bin_results
    
    # Tests
    
    def chi_sq_test(self, alpha):
    
        """
        Perform chi-sq test for different categories. 
    
        The first parameter is a dictionary of binned results.
        The second parameter denotes alpha, i.e. 1 - significance level for one-sided chi-sq test.
        """
    
        if alpha <0 or alpha>1:
            raise AssertionError("Alpha has to be from 0 to 1!")
    
        if not isinstance(self.bins,dict):
            raise AssertionError("Input has to be a Python dictionary!")
    
        ### This function assumes uniform bins being created.
        from scipy import stats
        exp_results = sum([self.bins[i] for i in self.bins])/len(self.bins)
    
        ### Calculate Chi-sq Statistic
        chi_sq =  sum([((self.bins[i] - exp_results)**2)/exp_results for i in self.bins])
    
        ### Set Significance Level and Critical Chi-Sq
        sig = 1 - alpha
        df = len(self.bins) - 1
        critical_chisq = stats.chi2.ppf(sig , df)
    
        ### Check for Significance Level
        chisq_test = chi_sq <= critical_chisq
    
        ### Prepare Statement
    
        if not chisq_test:
            stmt = '''
            Since the chi-sq value is larger than critical chi-sq value, the null hypothesis H0 
            that the random numbers generated are uniformly distributed is rejected. 
        
            \u03A7\N{SUPERSCRIPT TWO} = {0}
            Critical \u03A7\N{SUPERSCRIPT TWO} = {1}
            \u03A7\N{SUPERSCRIPT TWO} Result = {2}
        
            Please use another random number generator.
            '''.format(round(chi_sq,2), round(critical_chisq,2), 
            "Passed" if chisq_test else "Failed")
        else:
            stmt = '''
            Since the chi-sq value is smaller than critical chi-sq value, there is insufficient evidence
            that the random numbers generated are not uniformly distributed. Therefore,
            the null hypothesis H0 is not rejected.
            
            \u03A7\N{SUPERSCRIPT TWO} = {0}
            Critical \u03A7\N{SUPERSCRIPT TWO} = {1}
            \u03A7\N{SUPERSCRIPT TWO} Passed = {2}
        
            The random numbers passed the goodness-of-fit test.
            '''.format(round(chi_sq,2), round(critical_chisq,2), 
            "Passed" if chisq_test else "Failed")
    
        print(stmt)
        
    def up_and_down(self, alpha):
    
        from math import sqrt
        from scipy.stats import norm
    
        """
        The first parameter accepts a list of uniform random numbers. 
        The length of the numbers must be larger than 20.
    
        The second parameter denotes alpha, i.e. 1 - significance level for normal distribution.
        """
    
        if alpha <0 or alpha>1:
            raise AssertionError("Alpha has to be from 0 to 1!")
    
        if not isinstance(self.random_number,list):
            raise AssertionError("Input has to be a List!")
        
        delta = []
    
        for i in range(1, len(self.random_number)):
            if (self.random_number[i] - self.random_number[i-1])>0:
                delta.append(True)
                # True denotes increasing value.
            else:
                delta.append(False)
    
        up_and_down = []
        split_point = 0
    
        for i in range(1, len(delta)):
            if i == (len(delta)-1):
                up_and_down.append(delta[split_point:])
            elif delta[i-1] != delta[i]:
                up_and_down.append(delta[split_point:i])
                split_point = i
    
        A = len(up_and_down)
        EM = (2*len(self.random_number)-1)/3
        EV = (16*len(self.random_number)-29)/90
        Z0 = (A - EM)/sqrt(EV)
        Z_critical = norm.ppf(1-alpha/2)
        up_and_down_test = abs(Z0)<=Z_critical
    
        ### Prepare Statement
    
        if not up_and_down_test:
            stmt = '''
            Since the calculated Z value from the sample is larger than critical Z value, 
            the null hypothesis H0 that the random numbers generated are independent is rejected. 
        
            Z-Value = {0}
            Critical Z-Value = {1}
            Up-and-Down Test Result = {2}
        
            Please use another random number generator.
            '''.format(round(Z0,2), round(Z_critical,2), 
            "Passed" if up_and_down_test else "Failed")
        else:
            stmt = '''
            Since the calculated Z value from the sample is less than critical Z value, 
            the null hypothesis H0 that the random numbers generated are independent is not rejected. 
            There is insufficient evidence to show that the random numbers generated are not independent.
        
            Z-Value = {0}
            Critical Z-Value = {1}
            Up-and-Down Independence Test Result = {2}
        
            The random numbers passed the up-and-down independence test.
            '''.format(round(Z0,2), round(Z_critical,2), 
            "Passed" if up_and_down_test else "Failed")
    
        print(stmt)
    
    def above_and_below(self, alpha):
    
        from math import sqrt
        from scipy.stats import norm
    
        """
        The first parameter accepts a list of uniform random numbers.     
        The second parameter denotes alpha, i.e. 1 - significance level for normal distribution.
        """
    
        if alpha <0 or alpha>1:
            raise AssertionError("Alpha has to be from 0 to 1!")
    
        if not isinstance(self.random_number,list):
            raise AssertionError("Input has to be a List!")
    
        mean_test = [True if i>= 0.5 else False for i in self.random_number]
        above_and_below = []
        split_point = 0
    
        for i in range(1, len(mean_test)):
            if i == (len(mean_test)-1):
                above_and_below.append(mean_test[split_point:])
            elif mean_test[i-1] != mean_test[i]:
                above_and_below.append(mean_test[split_point:i])
                split_point = i
    
        B = len(above_and_below)
        n = len(mean_test)
        n1 = sum(mean_test)
        n2 = n - n1

        EM = 2*n1*n2/n + 0.5
        EV = (2*n1*n2*(2*n1*n2-n))/((n**2)*(n-1))
    
        Z0 = (B - EM)/sqrt(EV)
        Z_critical = norm.ppf(1-alpha/2)
        above_and_below_test = abs(Z0)<=Z_critical
    
        ### Prepare Statement
    
        if not above_and_below_test:
            stmt = '''
            Since the calculated Z value from the sample is larger than critical Z value, 
            the null hypothesis H0 that the random numbers generated are independent is rejected. 
        
            Z-Value = {0}
            Critical Z-Value = {1}
            "Above and Below the Mean" Test Result = {2}
        
            Please use another random number generator.
            '''.format(round(Z0,2), round(Z_critical,2), 
            "Passed" if above_and_below_test else "Failed")
        else:
            stmt = '''
            Since the calculated Z value from the sample is less than critical Z value, 
            the null hypothesis H0 that the random numbers generated are independent is not rejected. 
            There is insufficient evidence to show that the random numbers generated are not independent.
        
            Z-Value = {0}
            Critical Z-Value = {1}
            Above-and-Below the Mean Test Result = {2}
            
            The random numbers passed the "Above and Below the Mean" independence test.
            '''.format(round(Z0,2), round(Z_critical,2), 
            "Passed" if above_and_below_test else "Failed")
    
        print(stmt)

    def correlation_test(self, alpha):
    
        from math import sqrt
        from scipy.stats import norm
    
        """
        The first parameter accepts a list of uniform random numbers.     
        The second parameter denotes alpha, i.e. 1 - significance level for normal distribution.
        """
    
        if alpha <0 or alpha>1:
            raise AssertionError("Alpha has to be from 0 to 1!")
    
        if not isinstance(self.random_number,list):
            raise AssertionError("Input has to be a List!")
    
        rho_sum = 0
        n = len(self.random_number)
    
        for i in range(1,len(self.random_number)):
            number1 = self.random_number[i-1]
            number2 = self.random_number[i]
            rho_sum += number1*number2
    
        rho = 12/(n-1)*rho_sum-3
    
        EV = (13*n-19)/((n-1)**2)
    
        Z0 = rho/sqrt(EV)
        Z_critical = norm.ppf(1-alpha/2)
        corr_test = abs(Z0)<=Z_critical
    
        ### Prepare Statement
    
        if not corr_test:
            stmt = '''
            Since the calculated Z value from the sample is larger than critical Z value, 
            the null hypothesis H0 that the random numbers generated are independent is rejected. 
        
            Z-Value = {0}
            Critical Z-Value = {1}
            Correlation Test Result = {2}
        
            Please use another random number generator.
            '''.format(round(Z0,2), round(Z_critical,2), 
            "Passed" if corr_test else "Failed")
        else:
            stmt = '''
            Since the calculated Z value from the sample is less than critical Z value, 
            the null hypothesis H0 that the random numbers generated are independent is not rejected. 
            There is insufficient evidence to show that the random numbers generated are not independent.
        
            Z-Value = {0}
            Critical Z-Value = {1}
            Correlation Test Result = {2}
        
            The random numbers passed the Correlation Independence test.
            '''.format(round(Z0,2), round(Z_critical,2), 
            "Passed" if corr_test else "Failed")
    
        print(stmt)
        
    def plot_2d(self):
        if not isinstance(self.random_number,list):
            raise AssertionError("The input must be a list!")
    
        import matplotlib.pyplot as plt
    
        x = self.random_number
        y = x.copy()
        y.pop()
    
        plt.plot(x[1:], y, 'bo')
        plt.show()
    
    def plot_3d(self):
        if not isinstance(self.random_number,list):
            raise AssertionError("The input must be a list!")
    
        import matplotlib.pyplot as plt
    
        x = self.random_number
        y = x.copy()
        y.pop()
        z = y.copy()
        z.pop()
    
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax = plt.axes(projection='3d')
    
        ax.plot3D(x[2:], y[1:], z, 'bo')

        plt.show()

