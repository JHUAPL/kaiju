MPI Test Template
=================

.. code-block:: none

    module <INSERT MODULE NAME>
       use testHelperMpi
       # ADD ADDITIONAL USE STATEMENTS AS REQUIRED

       implicit none

       # DEFINE ANY OBJECTS THAT YOU WANT TO USE IN THE TESTS HERE
       # I RECOMMEND MAKING THEM ALLOCATABLE SO THAT THEY CAN BE DESTROYED BETWEEN TESTS

   contains

       @before
       subroutine setup(this)
           class (MpiTestMethod), intent(inout) :: this

           # THIS SHOULD ALLOCATE AND INITIALIZE ANY OBJECTS YOU WANT TO USE IN YOUR TESTS

           # YOU CAN GET THE MPI COMMUNICATOR CREATED FOR THIS TEST WITH THE FUNCTION
           # "getMpiF08Communicator(this)"

       end subroutine

       @after
       subroutine teardown(this)
           class (MpiTestMethod), intent(inout) :: this

           # THIS SHOULD DEALLOCATE AND DESTROY ANY OBJECTS USED IN YOUR TESTS
       end subroutine

       # IF YOU WANT TO CREATE ANY HELPER SUBROUTINES THAT YOUR TESTS USE, YOU CAN DO SO HERE

       @test(npes=[1,4,16,64])
       subroutine exampleTest(this)
           class (MpiTestMethod), intent(inout) :: this

           # THIS IS AN EXAMPLE TEST
           # THIS SHOULD BE DELETED, AND OTHER REAL TESTS SHOULD REPLACE IT BY PLACING THE
           # "@test" MARKER BEFORE THE TEST SUBROUTINES

           # UNLIKE THE SERIAL TEST METHODS, THE MPI METHODS ALSO HAVE THE "npes=([...])" SECTION
           # IN THE "@test" MARKER. THIS TELLS PFUNIT HOW MANY MPI RANKS TO RUN THE TEST WITH
           # IN THIS CASE, THIS EXAMPLE TEST WILL BE RUN 4 DIFFERENT TIMES, WITH
           # FIRST 1, THEN 4, THEN 16, AND THEN 64 MPI RANKS
           # THIS ALLOWS YOU TO CHECK DIFFERENT MPI SETUPS WITH THE SAME TEST

           # YOU CAN GET THE MPI COMMUNICATOR CREATED FOR THIS TEST WITH THE FUNCTION
           # "getMpiF08Communicator(this)"

       end subroutine

   end module
