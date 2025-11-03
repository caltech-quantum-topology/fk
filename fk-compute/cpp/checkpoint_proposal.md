
Build a checkpointing feature to start and stop the computation while it's running.

You do this by making a struct checkpointState and a class checkpointManager.
checkpointState will hold:
  an integer which counts how many times the computation has stopped and is resumed. This should be initialized to zero and incremented when a checkpointState is loaded by the checkpointManager
  an fkComputationEngine which gets saved every time a snapshot of the computation is taken
  a list of points, over which the computation needs to be carried through
  a the index of the last point that was compuited

checkpointManager should have two constructors. The first, should be declared with a checkpointState. It should create and store a unique name which will be used in the checkpoint file, and it should save the checkpointState immediately. The second, should be declared with a string representing the address of a saved checkpoint. It should load the checkpoint, save the name and it should return it. There should be a "save" method which takes a checkpointState and returns True if the saving was successfull or not.

Then, you should add some methods to fk_computation. On the user side, you should add a "computeWithCheckpointing" which takes a period of how often do you want to checkpoint. computeWithCheckpointing should do what compute does, but in addition it should also save a checkpoint every period. Then there should be a "resumeFromCheckpointing". This should take in a string (the location of the checkpoint file) and initialize the fk_computation from that state.

On the backend side fk_computation should have a saveCheckpoint method, where it saves the current state, calling checkpointManager::save, it should have a loadState which calls checkpointManger's loading constructor. Then it should have a "resume Computation" which resumes the computation.

Possibly the checkpoint should be saved as a binary so data structures can be loaded and saved with ease (think of pickle for python)
