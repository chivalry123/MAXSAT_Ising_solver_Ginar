#ifndef SERIALIZE_TUPLE
#define SERIALIZE_TUPLE

#include <boost/tuple/tuple.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/preprocessor/repetition.hpp>
#include <sstream>
#include <iostream>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>

namespace boost { namespace serialization {

    template<typename Archive, typename T1>

    void serialize(Archive & ar,
                   boost::tuples::cons<T1,boost::tuples::null_type> & t,

                   const unsigned int)
    {
      ar & boost::serialization::make_nvp("head",t.head);
    }

    template<typename Archive, typename T1, typename T2>

    void serialize(Archive & ar,
                   boost::tuples::cons<T1,T2> & t,

                   const unsigned int)
    {
      ar & boost::serialization::make_nvp("head",t.head);

      ar & boost::serialization::make_nvp("tail",t.tail);
    }

    template<typename Archive, typename T1>
    void serialize(Archive & ar,

                   boost::tuple<T1> & t,
                   const unsigned int)
    {
      ar & boost::serialization::make_nvp("head",t.head);
    }

#define GENERATE_TUPLE_SERIALIZE(z,nargs,unused)                            \
    template< typename Archive, BOOST_PP_ENUM_PARAMS(nargs,typename T) > \
    void serialize(Archive & ar,                                        \
                   boost::tuple< BOOST_PP_ENUM_PARAMS(nargs,T) > & t,   \
                   const unsigned int version)                          \
    {                                                                   \
      ar & boost::serialization::make_nvp("head",t.head);               \
      ar & boost::serialization::make_nvp("tail",t.tail);               \
    }


    BOOST_PP_REPEAT_FROM_TO(2,8,GENERATE_TUPLE_SERIALIZE,~);

}}

template<typename TupleType>

void test(TupleType t)
{
  std::ostringstream os;
  {
    boost::archive::text_oarchive ar(os);

    ar & boost::serialization::make_nvp("tuple",t);
  }
  TupleType t2;
  {

    std::istringstream is(os.str());
    boost::archive::text_iarchive ar(is);

    ar & boost::serialization::make_nvp("tuple",t2);
  }
  std::cout << t << " == " << t2  << " ? " << std::boolalpha << (t==t2) << std::endl;
}

#endif
