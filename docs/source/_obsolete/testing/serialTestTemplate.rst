Serial Test Template
====================

.. code-block:: bash

   module <INSERT MODULE NAME>
       use testHelper
       ! ADD ADDITIONAL USE STATEMENTS AS REQUIRED

       implicit none

       ! DEFINE ANY OBJECTS THAT YOU WANT TO USE IN THE TESTS HERE
       ! I RECOMMEND MAKING THEM ALLOCATABLE SO THAT THEY CAN BE DESTROYED BETWEEN TESTS

   contains

       @before
       subroutine setup()
           ! THIS SHOULD ALLOCATE AND INITIALIZE ANY OBJECTS YOU WANT TO USE IN YOUR TESTS
       end subroutine

       @after
       subroutine teardown()
           ! THIS SHOULD DEALLOCATE AND DESTROY ANY OBJECTS USED IN YOUR TESTS
       end subroutine

       ! IF YOU WANT TO CREATE ANY HELPER SUBROUTINES THAT YOUR TESTS USE, YOU CAN DO SO HERE

       @test
       subroutine exampleTest()
           ! THIS IS AN EXAMPLE TEST
           ! THIS SHOULD BE DELETED, AND OTHER REAL TESTS SHOULD REPLACE IT BY PLACING THE
           ! "@test" MARKER BEFORE THE TEST SUBROUTINES
       end subroutine

   end module
