!*********************************************************************
!**                                                                 **
!** "io_xml_input.F90"                                              **
!**                                                                 **
!** DESCRIPTION: Input class for XML formatted input. This class    **
!**              opens an XML formatted file and provides access to **
!**              parameters in the file                             **
!**                                                                 **
!** HISTORY                                                         **
!** Version  Date Changes                                Programmer **
!** 1.00 17/09/21 Original, ready for testing              J.L.Dean **
!** 1.10 17/09/28 Now uses kdefs                                JLD **
!** 1.20 17/10/04 Uses strings, verbose flag                    JLD **
!** 1.21 17/10/05 Get_Key_Val() is now a subroutine             JLD **
!** 1.22 18/01/29 Fixed Read() bug in Set_Str()                 JLD **
!** 1.23 18/02/28 Fixed empty line bug, log output is formatted JLD **
!** 1.30 18/03/14 Major updates to Is_Matched(),Get_Key_Val(),      **
!**               and Read_File()                               JLD **
!** 1.31 18/04/05 Added Exists()                                JLD **
!*********************************************************************

module XML_Data
   use kdefs
   use strings
   implicit none

   type XML_Data_T
      character(len=strLen) :: root ! Parent section of xml block
      character(len=strLen) :: sub ! sub-section containing key/value pairs
      character(len=strLen) :: pth ! the 'path' to the key/value pair. 
      character(len=strLen) :: key
      character(len=strLen) :: val
   end type XML_Data_T

   public :: Is_Equal

   ! Description:
   !   This function checks if a different XML_Data object is
   !   'equivalent' to self. Meaning if root, sub, key, and val
   !   of the XML_Data object have identical strings.
   !   Returns true if they are identical, False otherwise
   !
   ! Example Use:
   !   Check to see if xmld2 is identical to xmld1
   !
   ! ...
   ! type(XML_Data_T)  :: xmld1, xmld2
   ! logical         :: b
   ! ...
   ! b = kmld1%Is_Equal(xmld2)
   ! ...

   public :: Is_Matched

   ! This function is *very* similar to Is_Equal(),
   ! except there is no check on val. So, this returns
   ! true if root, sub, and key are the same, and false
   ! otherwise

contains

   logical function Is_Equal(this, xmld)
      class(XML_Data_T), intent(in) :: this
      type(XML_Data_T), intent(in)  :: xmld

      Is_Equal = Is_Matched(this, xmld)
      if (trim(toUpper(this%val)) /= trim((toUpper(xmld%val)))) then
         Is_Equal = .false.
      endif
   end function Is_Equal

   !-------------------------------------------------------------------

   logical function Is_Matched(this, xmld)
      class(XML_Data_T), intent(in) :: this
      type(XML_Data_T), intent(in)  :: xmld

      Is_Matched = .true.
      if (trim(toUpper(this%root)) /= trim(toUpper(xmld%root))) then
         Is_Matched = .false.
!      else if (trim(toUpper(this%sub)) /= trim(toUpper(xmld%sub))) then
!         Is_Matched = .false.
      else if (trim(toUpper(this%key)) /= trim(toUpper(xmld%key))) then
         Is_Matched = .false.
      endif
   end function Is_Matched

   !-------------------------------------------------------------------

end module XML_Data

!---------------------------------------------------------------------
!---------------------------------------------------------------------

module XML_Input
   use, intrinsic :: iso_fortran_env
   use XML_Data
   use kdefs

   implicit none

   logical, private :: doQuiet = .false.

   !Hard-coded format strings
   character(len=strLen), parameter, private :: realFormat = '(es12.5)'
   character(len=strLen), parameter, private :: strFormat  =  '(a12)'
   character(len=strLen), parameter, private :: boolFormat =  '(l12)'
   character(len=strLen), parameter, private :: intFormat  =  '(i12)'

   !Global strings, for default/file messages
   character(len=*), parameter :: defXML_Str = "Using DEFAULT value for "
   integer, parameter :: MaxXML = 1000 ! HARD-CODED!!

   private

   public :: XML_Input_T,SetAllXML2Quiet

   ! Description:
   !   The XML_Input_T "class". This "class" is designed to handle
   !   opening and reading of an XML formatted file. After successful
   !   read, the Get_* functions will be available to the User to
   !   extract data from the input file.
   !
   ! Example use of this class:
   !   ...
   !   use XML_Input
   !   type(XML_Input_T) :: inp
   !   real :: r,v
   !   logical :: verbosity
   !   integer :: p
   !   ...
   !   verbosity = .true.
   !   inp = New_XML_Input('input.XML','rootpath',verbosity)
   !   ...
   !   inp%Set_Val(r,'sect1/val1',0.05)
   !   inp%Set_Val(p,'/root2/sect1/val1',4)
   !   inp%Set_Val(v,'/root2/sect4/val1',24)
   !   ...
   !
   ! Note:
   !   See at Set_Val() for more information on getting/setting data
   !   and New_XML_Input() for how to initialize an object of this
   !   class.

   public :: New_XML_Input

   ! Description:
   !   Initialize the XML_Input 'object'. This function takes two
   !   arguments, the first being the file name (as string),
   !   the other being the root section of the xml file to read (also
   !   as string).
   !
   !
   ! Input Arguments:
   !   fname (Optional): the string representing the filepath and filename
   !                     of the xml file to be opened. If not provided,
   !                     "Input.XML" is set as the default filename
   !
   !   root (Optional): The default root section of the xml file (as string).
   !                    If not provided, no root section is defined. Read
   !                    "Set_Val()" for more information
   !
   !   verb (Optional): Verbosity flag for key/value pairs. If set to true,
   !                    all Set_* functions will write the key/value pairs
   !                    to terminal out. If not provided, default is false.
   !                    This is useful for logging/debugging purposes
   !
   !   For Example, given file 'Input.XML', and within it contains
   !   the following:
   !
   !   <root1>
   !     <sec>
   !       <subsec1 sect1 = "mySect" val1 = "val1" val2 = "val2"/>
   !       <subsec2 valA = "valA" valB = "valB" valC = "valC"/>
   !     </sec>
   !   </root1> 
   !   <root2>
   !     <sect1 = "mySect" val1 = "val1" val2 = "val2">
   !   </root2>
   !
   ! Initialize an XML_Input 'object' using the following example
   !   ...
   !   use class_XML_Input
   !   type(XML_Input_T) :: inp1
   !   ...
   !   inp1 = New_XML_Input('input.xml','root1')
   !   inp2 = New_XML_Input('myxml.xml','root2',.true.)
   !   ...
   !
   ! Now, Set_* functions can be called. (See Get_Val, for example)
   !
   ! Note:
   !   Both arguments are optional. If no arguments are
   !   provided, 'Input.XML' is set as the default file name
   !   and no 'root' path is defined.

   public :: ReadXmlImmediate

   !-------------------------------------------------------------------

   type XML_Input_T
      private
      character(len=strLen) :: mandatoryRoot ! required root of XML file and requests
      character(len=strLen) :: root ! default XML Root string
      character(len=strLen) :: fname ! fname
      character(len=strLen) :: buf ! string buffer
      type(XML_Data_T), dimension(:), allocatable :: xmld ! XML buffer
      integer :: fid = 0 ! ID number of opened File
      integer :: fst = 0 ! status of file
      logical :: vrb = .false.! flag for debug verbosity

      !-------------------------------------------------------------------


   contains
      procedure :: Open_File

      ! Private (internal) function. This function opens the file
      ! provided in New_XML_Input(), and reads/loads the XML content
      ! to an internal buffer. The user does not need to use this
      ! function explicitly as it is called internally by
      ! New_XML_Input() subroutine.
      !
      ! Returns True if the read was successful, False otherwise
      !
      ! Note:
      !   This is an internal function and should not be called
      !   by the user separately. However, this may become a
      !   public function in the future

      procedure :: Close_File

      ! Private (internal) function. This closes the file that is
      ! opened internally. The user does not need to use this
      ! function, and thus is set privately.
      !
      ! Returns True if file was closed successfully, False otherwise
      !
      ! Note:
      !   This is an internal function and should not be called
      !   by the user separately. However, this may become a
      !   public function in the future

      procedure :: Read_File

      ! Usage:
      !   if (inp%Read_File()) then ...
      !
      ! Description:
      !   This function performs the actual reading of the file;
      !   loading the xml contents into an internal buffer.
      !   Returns True if successful, false otherwise
      !
      ! Note:
      !   This is an internal function and should not be called
      !   by the user separately. However, this may become a
      !   public function in the future

      procedure :: Read_Line

      ! Usage:
      !   if (inp%Read_Line()) then ...
      !
      ! Description:
      !   Reads the next line of the input file into a buffer.
      !   Returns True if successful, false otherwise
      !
      ! Note:
      !   This is an internal function and should not (can not?)
      !   be called by the user. Returning false may indicate
      !   a problem with the input stream, or simply the end of
      !   input from that stream

      procedure :: Exists

      ! Usage:
      !   if (inp%Exists(path)) then ...
      !
      ! Description:
      !   This function checks the provided string for a match in any 
      !   of the XML_Data objects stored in self. Returns TRUE if found,
      !   FALSE otherwise.
      !
      ! Input Arguments:
      !   path:  char string representing the XML path to check for existence.
      !          str does *not* need to be the 'complete' path to a key/value
      !          pair. That is, partial paths are also checked.
      !
      ! Example Usage of this function:
      !   ...
      !   real :: v,w,x
      !   integer :: n,p
      !   ... (after initializing XML_Input inp) ...
      !   ...
      !   if (inp%Path_Exists("Remix/North")) then ... 
      !   if (inp%Path_Exists(str)) then ...
      !   if (inp%Path_Exists("Gamera/sim/runid")) then ... 
      !   ....
      !
      ! Note: 
      !   This function performs simple string comparisons for all 
      !   XML_Data objects stored in self. 
 
      procedure :: Get_Key_Val

      ! Usage:
      !   'Call inp%Get_Key_Val(str,res)'
      !
      ! Description:
      !   Searches for the matching key in the key/value pair within
      !   the internally stored xml_data variable, and returns the
      !   value as string. If there is no match, an empty
      !   string is returned
      !
      ! Input arguments:
      !   str: the string representing the xml path to the key of
      !        the key/value pair. As in, /Root/sect/key=value
      !
      !   res: The resulting string, containing the 'value' of the
      !        key/value pair. If the key path is not found, res
      !        will be an empty string.
      !
      ! Note:
      !   This function is called internally. The user does not need
      !   to use this function.

      !These map to RP (set in kdefs) from double or single precision
      procedure :: Set_Real_RpSp
      procedure :: Set_Real_RpDp
      ! procedure :: Set_Real_Sp

      ! ! See Set_Val. This is for 32-bit (single-precision?) real values

      ! procedure :: Set_Real_Dp

      ! ! See Set_Val. This is for 64-bit (double-precision?) real values

      procedure :: Set_Int

      ! See Set_Val.

      procedure :: Set_Bool

      ! See Set_Val.

      procedure :: Set_Str

      ! See Set_Val. Note that this will return the full string of the value in
      ! the key/value pair! 

      generic :: Set_Val => Set_Int, Set_Bool, Set_Str,Set_Real_RpDp,Set_Real_RpSp!Set_Real_Dp, Set_Real_Sp, 

      ! Usage:
      !   'call inp%Set_Val(v, str, def)'
      !
      ! Description:
      !   This sets v to the value associated with they key/value pair
      !   identified in str. If that key/value pair is not found, then
      !   v is set to def (being the default value)
      !
      ! Input arguments:
      !   v:   value to be assigned, the type is dictated by
      !        the default value type, given in def
      !
      !   str: the xml path to the value key. This str can take
      !        two forms; "/root/section/subsection/key" or "section/key".
      !        If the latter is used, a default string for 'root'
      !        which was provided within New_XML_Input(), is assumed
      !        otherwise it is not (ie, using the former example)
      !
      !   def: The default value to give should a value not be found
      !        for the xml path given in str. This also dictates
      !        the type of value of which v becomes. Thus, if
      !        def is of type int, v will also be int. (and so on)
      !
      ! Example Usage of this function:
      !   ...
      !   real :: v,w,x
      !   integer :: n,p
      !   ... (after initializing XML_Input inp) ...
      !   ...
      !   call inp%Set_Val(v,"/root1/sect2/val3",0.5)
      !   call inp%Set_Val(n,"sect3/val1",1)
      !   call inp%Set_Val(w,"/root2/sect2/val3",4)
      !   ....
      !
      ! Note:
      !   This function being generic allows for function
      !   overloading, meaning that the type of value returned
      !   is (in this case) dictated by the value type of the default
      !   value, as described above.
      !   The other Set_* functions (Set_Real_Dp, Set_Real_Sp, etc)
      !   are called depending on the value type of def. So, if def is
      !   a boolean, then internally, Set_Bool is called.
      !
      !   This subroutine should be the most widely used for users as
      !   it does not care about the type the user wishes to set.
      !   Errors will be thrown should a problem arise in dectecting the type
      procedure :: GetIsVerbose
      procedure :: SetVerbose
      procedure :: BeQuiet

      procedure :: GetFileStr
      procedure :: GetRootStr
      procedure :: SetRootStr

   end type XML_Input_T

contains

   subroutine GetFileStr(this,fStr)
      class(XML_Input_T), intent(in) :: this
      character(len=*), intent(inout)  :: fStr

      fStr = this%fname
   end subroutine GetFileStr

   subroutine GetRootStr(this,fStr)
      class(XML_Input_T), intent(in) :: this
      character(len=*), intent(inout)  :: fStr

      fStr = this%root
   end subroutine GetRootStr

   subroutine SetRootStr(this,fStr)
      class(XML_Input_T), intent(inout) :: this
      character(len=*), intent(in)   :: fStr

      this%root = trim(fStr)
   end subroutine SetRootStr

   subroutine GetIsVerbose(this, vrbOut)
      class(XML_Input_T), intent(in) :: this
      logical, intent(inout) :: vrbOut

      vrbOut = this%vrb
   end subroutine GetIsVerbose

   subroutine SetVerbose(this, setVrb)
      class(XML_Input_T), intent(inout) :: this
      logical, intent(in) :: setVrb

      this%vrb = setVrb
   end subroutine SetVerbose

   subroutine BeQuiet(this)
      class(XML_Input_T), intent(inout) :: this

      call SetVerbose(this, .false.)
   end subroutine BeQuiet

   subroutine Set_Real_RpSp(this, val, xmlp, dflt)
      class(XML_Input_T), intent(in) :: this
      real(rp), intent(inout)     :: val
      character(len=*), intent(in) :: xmlp
      real(sp), intent(in)        :: dflt
      character(len=strLen)       :: buf
      logical :: isXML

      val = dflt
      call this%Get_Key_Val(xmlp, buf)
      if (len(trim(buf)) /= 0) then
         isXML = .true.
         read (buf, *) val
      else
         isXML = .false.
      endif

      write(buf,realFormat) val
      call xmlOutput(this%root,xmlp,buf,this%vrb,isXML)

   end subroutine Set_Real_RpSp

   subroutine Set_Real_RpDp(this, val, xmlp, dflt)
      class(XML_Input_T), intent(in) :: this
      real(rp), intent(inout)     :: val
      character(len=*), intent(in) :: xmlp
      real(dp), intent(in)        :: dflt
      character(len=strLen)       :: buf
      logical :: isXML

      val = dflt
      call this%Get_Key_Val(xmlp, buf)
      if (len(trim(buf)) /= 0) then
         isXML = .true.
         read (buf, *) val
      else
         isXML = .false.
      endif

      write(buf,realFormat) val
      call xmlOutput(this%root,xmlp,buf,this%vrb,isXML)

   end subroutine Set_Real_RpDp


   subroutine Set_Int(this, val, xmlp, dflt)
      class(XML_Input_T), intent(in) :: this
      integer, intent(inout)       :: val
      character(len=*), intent(in) :: xmlp
      integer, intent(in)          :: dflt
      character(len=strLen)       :: buf
      logical :: isXML

      val = dflt
      call this%Get_Key_Val(xmlp, buf)
      if (len(trim(buf)) /= 0) then
         isXML = .true.
         read (buf, *) val
      else
         isXML = .false.
      endif
      write(buf,intFormat) val
      call xmlOutput(this%root,xmlp,buf,this%vrb,isXML)

   end subroutine Set_Int

   !------------------------------------------------------------------

   subroutine Set_Bool(this, val, xmlp, dflt)
      class(XML_Input_T), intent(in) :: this
      logical, intent(inout)       :: val
      character(len=*), intent(in) :: xmlp
      logical, intent(in)          :: dflt
      character(len=strLen)       :: buf
      logical :: isXML

      val = dflt
      call this%Get_Key_Val(xmlp, buf)
      if (len(trim(buf)) /= 0) then
         isXML = .true.
         read (buf, *) val
      else
         isXML = .false.
      endif

      write(buf,boolFormat) val
      call xmlOutput(this%root,xmlp,buf,this%vrb,isXML)

   end subroutine Set_Bool

   !------------------------------------------------------------------

   subroutine Set_Str(this, val, xmlp, dflt)
      class(XML_Input_T), intent(in) :: this
      character(len=*), intent(inout) :: val
      character(len=*), intent(in)    :: xmlp
      character(len=*), intent(in)    :: dflt
      character(len=strLen)          :: buf
      logical :: isXML

      val = dflt
      call this%Get_Key_Val(xmlp, buf)
      if (len(trim(buf)) /= 0) then
         isXML = .true.
         read (buf, *) val
      else
         isXML = .false.
      endif
      write(buf,strFormat) val
      call xmlOutput(this%root,xmlp,buf,this%vrb,isXML)

   end subroutine Set_Str

   !------------------------------------------------------------------

   subroutine xmlOutput(rStr,xStr,vStr,doVerb,isXML)
      character(len=*), intent(in) :: rStr,xStr,vStr
      logical, intent(in) :: doVerb,isXML
      character(len=strLen) :: bStr,srcStr
      character(len=*), parameter :: formt  = '(2X, a40, a1, 2X, a12, a10)'
      character(len=*), parameter :: formtC = '(2X, a40, a1, 2X, a12, a,a10,a)'
      integer :: pos

      if ( (.not. doVerb) .or. (doQuiet) ) return
      pos = scan(trim(xStr), '/')
      if (pos /= 1) then
          ! key is relative to root
          bStr = trim(toUpper(rStr)) // '/' // trim(xStr)
      else
          ! key is absolute
          bStr = trim(xStr)
      endif
      
      if (isXML) then
         srcStr = "(XML)"
         write(*,formtC) adjustl(bStr), ':', trim(vStr), ANSIBLUE, trim(srcStr), ANSIRESET
      else
         srcStr = "(DEFAULT)"
         write(*,formtC) adjustl(bStr), ':', trim(vStr), ANSIGREEN, trim(srcStr), ANSIRESET
      endif
      !write(*,formt) adjustl(bStr), ':', trim(vStr), trim(srcStr)
   end subroutine xmlOutput

   function New_XML_Input(fname, root, verbOpt, fileRoot)
      type(XML_Input_T) :: New_XML_Input
      character(len=*), intent(in) :: fname
      character(len=*), intent(in) :: root
      logical, intent(in), optional :: verbOpt
      character(len=*), intent(in), optional :: fileRoot

      integer :: fid
      character(len=strLen) :: buf

      fid = 123 ! - HARDCODED FILE ID!
      !   this parameter may need to be pushed
      !   onto the user such that they protect
      !   against file ID collisions!

      New_XML_Input%fid = fid
      New_XML_Input%root = root
      New_XML_Input%fname = fname

      New_XML_Input%vrb = .false.
      if (present(verbOpt)) then
         New_XML_Input%vrb = verbOpt
      endif

      if (present(fileRoot)) then
          New_XML_Input%mandatoryRoot = trim(fileRoot)
      else
          ! default mandatory root
          New_XML_Input%mandatoryRoot = "kaiju"
      endif

      ! Warn if provided an impty file/root string
      if (New_XML_Input%fname == "") then
         write (*, *) "WARNING: No filename provided."
         write (*, *) "         defaulting to Input.XML"
         New_XML_Input%fname = "Input.XML"
      endif
      if (root == "") then
         write (*, *) "WARNING: No XML Root provided."
      endif

      ! Throw an error if unable to Open the file
      if (.NOT. (New_XML_Input%Open_File())) then
         write (*, *) "!!!ERROR!!! Unable to Open File", New_XML_Input%fname
         stop
      endif


      ! Updated to new Read_File function -- JLD 18/03/09
      if (.NOT. (New_XML_Input%Read_File())) then
        write (*, *) "!!!ERROR!!! Error encountered while reading input file"
        stop
      endif

      ! Throw an error if unable to Close the file
      if (.NOT. (New_XML_Input%Close_File())) then
         write (*, *) "!!!ERROR!!! Unable to Close File"
         write (*, *) New_XML_Input%fname
         stop
      endif
   end function New_XML_Input

   !------------------------------------------------------------------

   logical function Open_File(this)
      class(XML_Input_T), intent(inout) :: this

      open (this%fid, file=trim(this%fname), status='old', action='read', iostat=this%fst)
      if (this%fst /= 0) then
         Open_File = .false.
      else
         Open_File = .true.
      end if
   end function Open_File

   !------------------------------------------------------------------

   logical function Close_File(this)
      class(XML_Input_T), intent(inout) :: this

      close (this%fid, iostat=this%fst)
      if (this%fst /= 0) then
         Close_File = .false.
      else
         Close_File = .true.
      endif
   end function Close_File

   !------------------------------------------------------------------

   logical function Read_File(this)
     class(XML_Input_T), intent(inout) :: this
     integer               :: pos
     integer               :: nosg
     integer               :: nxml
     character(len=strLen) :: pth,rt,sec,key,val,tmp
     type(XML_Data_T), dimension(MaxXML) :: xmld

     nxml = 0
     nosg = 0 ! number of 'open' segments

     do while (this%Read_Line())

       ! skip this loop if current line is an empty line
       if(len_trim(this%buf)==0) then
         cycle
       endif

       ! Check (and skip) Comments; starting with <! or <?
       if(index(trim(this%buf), '<!') /= 0 .or. &
       &  index(trim(this%buf), '<?') /= 0) then 
         cycle
       endif

       ! Possible combinations per line: 
       !   "< ... >"         <--- Root or Section
       !   "</ ... >"        <--- End of root or Section
       !   "< ... = ... />"  <--- section with key/val pairs
       !   "< ... = ... >"   <--- same as previous

       ! First, check for the following; "<", and *not* "</" or "="
       ! This indicates either a new root or a new section
       if (index(trim(this%buf), '<')  /= 0 .and. &
       &   index(trim(this%buf), '</') == 0 .and. &
       &   index(trim(this%buf), '=')  == 0) then
       
         nosg = nosg + 1  ! Increment no. of 'open' sections
         
         ! Remove < and > from buffer
         this%buf = adjustl(trim(this%buf))
         pos = scan(this%buf, '<')
         this%buf = trim(this%buf(pos + 1:))
         pos = scan(this%buf, '>')
         this%buf = trim(this%buf(1:pos - 1))
       
         ! If nosg == 1, we have the root name:
         if(nosg == 1) then
           rt = this%buf  ! set root       
           pth = ""       ! and initialize the key 'path' string
 
         ! Otherwise; this is the opening of a new subsection
         else 
           ! First check if the 'pth' string is empty, 
           if(len(trim(pth)) .le. 0) then
             pth = trim(this%buf)

           ! If it is 'not' empty, append '/' and the next 
           ! section to pth 
           else 
             pth = trim(pth) // '/'
             pth = trim(pth) // trim(this%buf)
           endif
         endif

       ! If '</' is found; decrement the no. of 'open' sections, 
       ! and remove the 
       else if (index(trim(this%buf), '</') /= 0) then
         nosg = nosg - 1
         pos = scan(trim(pth),'/',.true.) ! Find the last '/' character

         if(pos > 0) then
           pth = trim(pth(1:pos-1))
         else
           pth = "" ! no longer in anything but root
         endif

       ! Finally, if '=' is found; it's a section containing key/value pairs.
       ! the 'key' will be, essentially 'sec/subsec/subsubsec/key'.  
       else if (index(trim(this%buf), '=') /= 0) then

         ! Remove the first '<'. 
         pos = scan(trim(this%buf),'<')
         this%buf = trim(this%buf(pos+1:))
          
         !The immediate next set of chars, (which must be followed by a space)
         !is the 'last' subsection of this path.   
         pos = scan(trim(this%buf), ' ')
         
         !Set a temporary container for this last subsection
         tmp = trim(this%buf(:pos-1))
        
         !And modify buffer such that only 'key="val"' pairs are left (plus a />
         !or >)

         this%buf = trim(this%buf(pos+1:))

         ! Now parse the rest of the sting searching for
         ! key/value pairs
         pos = 1 ! Initialize position value
         do

           pos = scan(this%buf, '=')
           if (pos == 0) then
              exit
           endif

           nxml = nxml+1
           if (nxml > size(xmld)-1) then
             write (*, *) "!!!ERROR!!! Total no. of input values exceeds temp xml data array size"
             stop
           endif

           ! First: set root
           xmld(nxml)%root = trim(rt)

           ! if nosg > 1, append contents of tmp into pth, and set
           ! tmp = path
           if (nosg > 1) then

             xmld(nxml)%key = trim(pth)
             xmld(nxml)%key = trim(xmld(nxml)%key) // '/'
             xmld(nxml)%key = trim(xmld(nxml)%key) // trim(tmp)
             xmld(nxml)%key = trim(xmld(nxml)%key) // '/'
             xmld(nxml)%key = trim(xmld(nxml)%key) // trim(this%buf(1:pos-1))
    
           ! Otherwise; append key to tmp. 
           else if (nosg == 1) then
             xmld(nxml)%key = trim(tmp)
             xmld(nxml)%key = trim(xmld(nxml)%key) // '/'
             xmld(nxml)%key = trim(xmld(nxml)%key) // trim(this%buf(1:pos-1))
           endif

           ! Don't set the sub -- just set key to be the 'path'
           ! variable tmp contains the sub portion 
           this%buf = trim(this%buf(pos+1:))
           pos = scan(this%buf, '"') ! remove first "
           this%buf = trim(this%buf(pos + 1:))
           pos = scan(this%buf, '"') ! get pos of second "
           xmld(nxml)%val = this%buf(1:pos - 1) ! set value
           this%buf = trim(this%buf(pos + 2:)) ! remove value"_

         end do

       endif
     end do

     ! Finally, check for any problems. 
     ! First problem would be if any sections are stil 'open' 
     ! (ie, if nosg > 0 )
     if (nosg > 0) then
       Read_File = .false.
     
     ! Otherwise! Allocate self and return true
     else 
       allocate (this%xmld, source=xmld(1:nxml))
       Read_File = .true.
     endif

  end function Read_File

   !------------------------------------------------------------------

   subroutine Get_Key_Val(this, str, sout, checkRoot)
      class(XML_Input_T), intent(in)   :: this
      character(len=*), intent(in)   :: str
      character(len=strLen), intent(inout) :: sout
      logical, intent(in), optional :: checkRoot

      integer            :: i
      integer            :: pos
      character(len=strLen) :: buf
      type(XML_Data_T)     :: xmld
      logical            :: doRootCheck

      if(present(checkRoot)) then
          doRootCheck = checkRoot
      else
          doRootCheck = .true. ! default
      endif

      sout = "" ! initialize the result

      ! check if string begins with '/' character. If it's not the
      ! first position, assume the given root path.
      buf = trim(str)
      pos = scan(trim(buf), '/')
      if (pos /= 1) then
         ! Check if the default root path is provided. If not, stop
         ! and return 0 (not found)
         if (this%root == "") then
            return
         endif
         xmld%root = this%root ! otherwise, set the root using default
         ! deal with the edge case that the root has more than entry in it
         pos = scan(trim(this%root), '/')
         if(pos .gt. 0) then
             ! the root contains a '/'
             xmld%root = this%root(1:pos-1);
             ! combine the rest of the root into the buffer that goes into the key
             buf = trim(this%root(pos+1:)) // '/' // trim(buf(:))
         endif
      else
         ! If however, it 'is' in the first position, parse
         ! the string and set root
         buf = trim(buf(pos + 1:))
         pos = scan(trim(buf), '/')
         xmld%root = buf(1:pos - 1)
         buf = trim(buf(pos + 1:))
      endif

      ! and finally, set the key string. This time, no '/'
      ! should be found, and trim(buf) should just be the key
      xmld%key = trim(buf)
    
      ! check if the search root matches the required root
      if(doRootcheck .and. len(this%mandatoryRoot) > 0) then
          if(ToUpper(xmld%root) /= ToUpper(this%mandatoryRoot)) then
              write(*,*) "ERROR:"
              write(*,*) "The root of requested XML keys MUST be '", trim(this%mandatoryRoot), "'"
              write(*,*) "The following key was requested with an incorrect root:"
              write(*,*) "  /",trim(xmld%root),"/",trim(xmld%key)
              write(*,*) "Stopping application."
              write(*,*) ""
              stop
          endif
      endif

      ! Now, search through the internal xml_data type
      ! to extract the value. If the value is not found,
      ! return an empty string
      do i = 1, size(this%xmld), 1
         if (Is_Matched(this%xmld(i), xmld)) then
            sout = this%xmld(i)%val
            return
         endif
      end do

   end subroutine Get_Key_Val

   !------------------------------------------------------------------

   logical function Read_Line(this)
      class(XML_Input_T), intent(inout) :: this
      integer          :: bsz
      character(len=strLen) :: buf

      read(this%fid, '(A)', iostat=this%fst) buf
      this%buf = trim(buf)
      if (this%fst /= 0) then
         Read_Line = .false.
      else
         Read_Line = .true.
      end if
   end function Read_Line

   !------------------------------------------------------------------

   logical function Exists(this, path)
      class(XML_Input_T), intent(in)  :: this
      character(len=*), intent(in)    :: path
      integer                         :: i, matchI
      character(len=strLen)           :: buf
      Exists = .false.
      if(size(this%xmld).eq.0) return
      do i=1, size(this%xmld), 1
         buf = trim(this%xmld(i)%root) // &
               trim("/") // trim(this%xmld(i)%key)
         matchI = index(toUpper(trim(buf)),toUpper(trim(path)))
         if(matchI .ge. 1 .and. matchI == (1 + len_trim(buf) - len_trim(path))) then
            ! allow the path to start inside if the user only gave part of the path
            ! but insist that it match the end of the buffer, so that we match the entire
            ! variable name
            Exists = .true.
            return
         endif
      end do
   end function Exists

   !------------------------------------------------------------------
   
   !Set all XML reading to quiet on this rank
   subroutine SetAllXML2Quiet()
      doQuiet = .true.
   end subroutine SetAllXML2Quiet

   !------------------------------------------------------------------

   !Helper function to open an xml, read a single parameter, and close it
   !For when you just need that one parameter, right now
   subroutine ReadXmlImmediate(fname, keyIn, valOut, dflt, verbOpt)
      character(len=*), intent(in) :: keyIn
      character(len=strLen), intent(inout) :: valOut
      character(len=*), intent(in) :: fname
      character(len=*), intent(in) :: dflt
      logical, intent(in), optional  :: verbOpt

      type(XML_Input_T) :: xmlInp

      if (present(verbOpt)) then
         xmlInp = New_XML_Input(trim(fname),'/',verbOpt)
      else
         xmlInp = New_XML_Input(trim(fname),'/')
      endif

      call xmlInp%Get_Key_Val(trim(keyIn), valOut)
      if(len(trim(valOut)) == 0) valOut = trim(dflt)

   end subroutine

end module XML_Input

!--------------------------------------------------------------------
