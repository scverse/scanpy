External tools
==============

For external tools we’d prefer the code mostly lives in an external package.
For internal, it’ll be a more in-depth review as we'd be more responsible for maintenance and would essentially be endorsing the method.

Testing
-------

Tests for tools in the external module are less thorough than tools in the main API.

While tests for the main API are there to assure correctness, we assume external tools have their own test suites in their package.
Instead, tests here are to let us know:

* If changes in scanpy break the external tool
* If changes in dependencies would break external tools
