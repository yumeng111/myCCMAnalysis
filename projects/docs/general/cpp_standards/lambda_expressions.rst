.. highlight:: cpp

Lambda Expressions
------------------

* "A lambda expression is just that: an expression. It is part of the source code" [Meyer]_.
* In other words, we can say that a lambda expression allows to write a function inside of another function.

``auto`` can be used as a ``lambda expression``, as a ``type deduction``, or to handle STL containers.

``auto`` as a lambda expression
"""""""""""""""""""""""""""""""

If you place the below expression inside of a ``main()`` function (for example), this will work::

  const auto lambda_momentum = [] (const double &px, const double &py){
    const double pz = px * 0.6 + py * 0.8;
    const double &momentum = std::sqrt(px * px + py * py + pz * pz);
    return momentum;
  };


``auto`` as type deduction
""""""""""""""""""""""""""

``auto`` can be used as a type. This is very useful when objects are
very large and hard to handle. See the ``auto`` documentation for more
explanation about type deduction.  For example:
``std::map<std::string,int>::iterator``::

  for (auto it=StringIntMap.begin(); it != StringIntMap.end(); it++){
        std::cout << "Key: " << it->first << " Value " << it->second << std::endl;
  }

In the above example, ``auto`` does not work as lambda expression, this works as type deduction.


Use a lambda expression to handle a STL container
"""""""""""""""""""""""""""""""""""""""""""""""""

You always should prefer to use STL algorithms instead of handwritten
calculations [Sutter_Alexandrescu]_. Below an example of how to use
lambda expressions to manipulate a STL vector::

  std::vector<int>::iterator first_even = // this line can be replaced by "auto", used as Type Deduction
    std::find_if(IntVect.begin(), IntVect.end(), [](int it) {
      return it % 2==0;
    });

  std::cout << "The first even number is " << * first_even << std::endl;

References
""""""""""

For further, more complete examples, please refer to `this code`_ from the August 2016 code sprint.

.. [Meyer] *Effective Modern C++* by Scott Meyers (O’Reilly). Copyright 2015 Scott Meyers, 978-1-491-90399-5.

.. [Sutter_Alexandrescu] *C++ coding standards : 101 rules, guidelines, and best practices* by Herb Sutter, Andrei Alexandrescu (Addison-Wesley Professional). Copyright 2005 Pearson Education, Inc, 0-321-11358-6.

.. _this code: http://code.icecube.wisc.edu/projects/icecube/browser/IceCube/sandbox/CodeSprints/August2016/lambdas.cc
